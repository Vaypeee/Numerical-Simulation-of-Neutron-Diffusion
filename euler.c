#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "euler.h"
#include "geometry.h"
#include "prob.h"
#include "residual.h"
#include "critical.h"
#include "mytime.h"
#include "interface_primme.h"

#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MKDIR(p) mkdir(p, 0755)
#endif

/* Facteurs d'echelle etudies, relatifs a la dimension critique. */
static const double ECHELLES[]  = { 0.98, 1.00, 1.02 };
static const int     N_ECHELLES = 3;

/* Fraction de la borne de stabilite utilisee comme pas de temps. */
#define SECURITE_DT 0.9

/* Nombre maximal de points ecrits dans le fichier de sortie par courbe. */
#define MAX_POINTS 2000

/* Nombre de pas de la demonstration de stabilite. Il en faut assez pour
   que meme une amplification par pas tres legerement superieure a 1 (par
   exemple 1.002 a dt = 1.001 dt_max) devienne visible, sachant que le
   mode concerne part avec une amplitude tres faible dans le flux initial. */
#define N_PAS_DEMO 20000

/*
   Borne de Gershgorin sur la plus grande valeur propre de A.
   ==========================================================
   Pour une matrice symetrique, toute valeur propre est contenue dans la
   reunion des disques de Gershgorin ; la borne superieure est donc

        lambda_max  <=  max_i ( a_ii + somme_{j != i} |a_ij| )

   c'est-a-dire le maximum des sommes de valeurs absolues par ligne. Pour le
   laplacien a cinq points cette borne vaut 8/h^2, atteinte asymptotiquement.
   On la calcule ici directement depuis la structure CSR, ce qui reste
   valable si la matrice change.
*/
static double gershgorin_max(int n, const int *ia, const double *a)
{
    double borne = 0.0;
    int    i, k;

    for (i = 0; i < n; i++) {
        double somme = 0.0;
        for (k = ia[i]; k < ia[i + 1]; k++)
            somme += fabs(a[k]);
        if (somme > borne)
            borne = somme;
    }

    return borne;
}

/*
   Un pas d'Euler progressif, fusionne.
   ====================================
        Phi^{k+1} = Phi^k - dt * tau * ( A Phi^k - beta^2 Phi^k )

   Comme pour le calcul du residu (tache 2), la ligne de A n'est parcourue
   qu'une fois et le resultat intermediaire ne quitte pas les registres : on
   n'ecrit jamais A*Phi dans un tableau temporaire. La norme du nouveau
   vecteur est accumulee dans la meme passe, ce qui evite une relecture
   complete de Phi^{k+1}.

   Retourne || Phi^{k+1} ||_2 et ecrit le maximum dans *vmax.
*/
static double euler_step(int n, const int *ia, const int *ja, const double *a,
                         double coef, double beta2,
                         const double *phi, double *phi_new, double *vmax)
{
    const int    * const restrict pia = ia;
    const int    * const restrict pja = ja;
    const double * const restrict pa  = a;
    const double * const restrict px  = phi;
    double       * const restrict py  = phi_new;

    double norme2 = 0.0, maxi = 0.0;
    int    i, k, k_deb = pia[0];

    for (i = 0; i < n; i++) {
        const int    k_fin = pia[i + 1];
        const double xi    = px[i];
        double       r     = -beta2 * xi;   /* (A - beta^2 I) Phi, ligne i */
        double       v;

        for (k = k_deb; k < k_fin; k++)
            r += pa[k] * px[pja[k]];

        v      = xi - coef * r;             /* coef = dt * tau */
        py[i]  = v;
        norme2 += v * v;
        if (v > maxi) maxi = v;

        k_deb = k_fin;
    }

    *vmax = maxi;
    return sqrt(norme2);
}

/*
   Integre le probleme de Cauchy pour un facteur d'echelle donne.
   ==============================================================
   Le reacteur de cote s * L_crit est obtenu sans remailler : d'apres la
   tache 3, dilater le domaine d'un facteur alpha divise la matrice par
   alpha^2. On multiplie donc simplement les coefficients de A par
   (L / L_reacteur)^2.

   Consequence utile : pour s = 1 exactement, la plus petite valeur propre
   de la matrice mise a l'echelle vaut lambda_min / alpha^2 = beta^2 a la
   precision machine. Le cas critique est donc EXACTEMENT critique au niveau
   discret, et le plateau observe n'est pas pollue par l'erreur de
   discretisation.

   Les valeurs de || Phi(t) || sont ecrites dans le fichier deja ouvert 'f',
   sous forme d'un bloc separe des autres par deux lignes blanches (index
   gnuplot).
*/
static int integre(FILE *f, int n, const int *ia, const int *ja,
                   const double *a_ref, const double *v1,
                   double lambda_ref, double alpha, double echelle,
                   double t_final, int index_bloc)
{
    double *a      = malloc((size_t) ia[n] * sizeof(double));
    double *phi    = malloc((size_t) n * sizeof(double));
    double *phi_np = malloc((size_t) n * sizeof(double));
    double  facteur, lambda_min, omega, dt, dt_max, lam_max;
    double  norme0, norme, vmax, t, cos_angle, ps, nv1;
    double  norme_demi, t_demi, omega_obs;
    long    pas, n_pas, periode;
    int     i, k;

    if (!a || !phi || !phi_np) {
        printf("\n ERREUR : pas assez de memoire pour l'integration\n\n");
        free(a); free(phi); free(phi_np);
        return 1;
    }

    /* --- matrice du reacteur de cote (echelle * L_crit) ------------------ */
    facteur = 1.0 / ((alpha * echelle) * (alpha * echelle));
    for (k = 0; k < ia[n]; k++)
        a[k] = a_ref[k] * facteur;

    lambda_min = lambda_ref * facteur;
    lam_max    = gershgorin_max(n, ia, a);

    /* --- taux de croissance theorique ------------------------------------
       Phi(t) = somme_i c_i exp(-tau (lambda_i - beta^2) t) v_i ; aux grands
       temps le mode fondamental domine et || Phi || ~ exp(omega t) avec  */
    omega = REACTOR.tau * (REACTOR.beta2 - lambda_min);

    /* --- pas de temps ----------------------------------------------------
       Euler progressif applique a dPhi/dt = -tau (A - beta^2 I) Phi est
       stable si et seulement si | 1 - dt tau (lambda_i - beta^2) | <= 1 pour
       toute valeur propre, soit dt <= 2 / (tau (lambda_max - beta^2)). */
    dt_max = 2.0 / (REACTOR.tau * (lam_max - REACTOR.beta2));
    dt     = SECURITE_DT * dt_max;
    n_pas  = (long) (t_final / dt) + 1;

    /* --- condition initiale : flux constant phi0 ------------------------- */
    for (i = 0; i < n; i++)
        phi[i] = REACTOR.phi0;

    norme0 = 0.0;
    for (i = 0; i < n; i++)
        norme0 += phi[i] * phi[i];
    norme0 = sqrt(norme0);

    printf("\n  echelle s = %.2f  ->  cote = %.6f m\n", echelle,
           echelle * alpha * REACTOR.L);
    printf("    lambda_min = %.10f m^-2   (beta^2 = %g)\n", lambda_min, REACTOR.beta2);
    printf("    regime : %s\n",
           fabs(omega) < 1e-9 ? "CRITIQUE (stationnaire)"
                              : (omega > 0.0 ? "SUR-CRITIQUE (croissance)"
                                             : "SOUS-CRITIQUE (decroissance)"));
    printf("    taux theorique omega = tau (beta^2 - lambda_min) = %+.6f s^-1\n", omega);
    printf("    dt_max (stabilite) = %.4e s   dt utilise = %.4e s   (%ld pas)\n",
           dt_max, dt, n_pas);

    /* --- boucle en temps -------------------------------------------------- */
    periode = n_pas / MAX_POINTS;
    if (periode < 1) periode = 1;

    fprintf(f, "# bloc %d : echelle s = %.2f, omega theorique = %.6f s^-1\n",
            index_bloc, echelle, omega);
    fprintf(f, "# t [s]   ||Phi(t)|| / ||Phi(0)||   max Phi [m^-2 s^-1]\n");
    fprintf(f, "0 1 %.8g\n", REACTOR.phi0);

    norme = norme0;
    norme_demi = norme0;
    for (pas = 1; pas <= n_pas; pas++) {
        norme = euler_step(n, ia, ja, a, dt * REACTOR.tau, REACTOR.beta2,
                           phi, phi_np, &vmax);

        /* echange des tampons : pas de recopie */
        { double *tmp = phi; phi = phi_np; phi_np = tmp; }

        /* norme a mi-parcours : sert a mesurer le taux de croissance sur la
           seconde moitie de la simulation, ou les modes eleves ont deja
           disparu et ou seul le mode fondamental subsiste. */
        if (pas == n_pas / 2)
            norme_demi = norme;

        if (pas % periode == 0 || pas == n_pas) {
            t = pas * dt;
            fprintf(f, "%.8g %.8g %.8g\n", t, norme / norme0, vmax);
        }

        if (!(norme > 0.0) || norme > 1e300) {   /* divergence ou NaN */
            printf("    ATTENTION : divergence detectee au pas %ld\n", pas);
            break;
        }
    }
    fprintf(f, "\n\n");   /* separateur de bloc pour gnuplot */

    /* --- le profil final est-il le mode fondamental ? --------------------
       On compare la direction de Phi(t_final) a celle du vecteur propre
       fondamental v1 : le cosinus de l'angle doit tendre vers 1, puisque
       tous les autres modes decroissent plus vite. */
    ps = 0.0; nv1 = 0.0;
    for (i = 0; i < n; i++) {
        ps  += phi[i] * v1[i];
        nv1 += v1[i] * v1[i];
    }
    cos_angle = fabs(ps) / (norme * sqrt(nv1));

    /* --- taux de croissance mesure --------------------------------------
       Sur la seconde moitie de la simulation, || Phi || ~ C exp(omega t) :
       le taux se lit donc directement dans le rapport des normes. */
    t_demi    = (n_pas / 2) * dt;
    omega_obs = log(norme / norme_demi) / (n_pas * dt - t_demi);

    printf("    amplification || Phi(t_f) || / || Phi(0) || = %.6e\n", norme / norme0);
    printf("    taux mesure  omega_obs = %+.6f s^-1   (theorique %+.6f, ecart %.2e)\n",
           omega_obs, omega, fabs(omega_obs - omega));
    /* Euler progressif amplifie chaque mode par (1 + dt*omega) au lieu de
       exp(dt*omega). Le taux mesure vaut donc log(1+dt*omega)/dt, soit
       omega - dt*omega^2/2 au premier ordre : l'ecart ci-dessus doit
       coincider avec cette erreur de troncature d'ordre 1. */
    printf("      erreur de troncature attendue (ordre 1) : dt*omega^2/2 = %.2e\n",
           0.5 * dt * omega * omega);
    printf("    projection du flux initial sur le mode fondamental :\n");
    printf("      || Phi(t_f) || / (|| Phi(0) || exp(omega t_f)) = %.6f\n",
           (norme / norme0) / exp(omega * t_final));
    printf("      (cette constante ne depend pas de s : c'est le poids du mode\n");
    printf("       fondamental dans le flux initial constant phi0)\n");
    printf("    cos(angle entre Phi(t_f) et le mode fondamental) = %.10f\n", cos_angle);

    free(a); free(phi); free(phi_np);
    return 0;
}

/*
   Ecrit le script gnuplot de la figure d'evolution.
*/
static int write_script(const char *nom_gp, const char *nom_dat,
                        const char *png, int m, double t_final)
{
    FILE *f = fopen(nom_gp, "w");
    int   k;

    if (f == NULL) {
        printf("\n ERREUR : impossible d'ecrire '%s'\n\n", nom_gp);
        return 1;
    }

    fprintf(f, "# Script genere automatiquement par euler.c\n");
    fprintf(f, "set terminal pngcairo size 1000,720 enhanced font 'Arial,13'\n");
    fprintf(f, "set output '%s'\n", png);
    fprintf(f, "set title \"Evolution du flux par Euler progressif"
               " (m = %d)\\ndPhi/dt = -{/Symbol t} (A - {/Symbol b}^2 I) Phi,"
               "  Phi(0) = {/Symbol f}_0 = %g m^{-2} s^{-1}\"\n", m, REACTOR.phi0);
    fprintf(f, "set xlabel \"temps t  [s]\"\n");
    fprintf(f, "set ylabel \"|| {/Symbol F}(t) || / || {/Symbol F}(0) ||\"\n");
    fprintf(f, "set logscale y\n");
    fprintf(f, "set grid\n");
    fprintf(f, "set key left top\n");
    fprintf(f, "set xrange [0:%.8g]\n", t_final);

    fprintf(f, "plot ");
    for (k = 0; k < N_ECHELLES; k++) {
        const char *couleur = (k == 0) ? "#2171b5" : (k == 1) ? "#000000" : "#a50026";
        const char *regime  = (ECHELLES[k] < 1.0) ? "sous-critique"
                            : (ECHELLES[k] > 1.0) ? "sur-critique" : "critique";
        fprintf(f, "%s'%s' index %d using 1:2 with lines lw 3 lc rgb '%s' "
                   "title 's = %.2f  (%s)'",
                k ? ", " : "", nom_dat, k, couleur, ECHELLES[k], regime);
    }
    fprintf(f, "\n");

    fclose(f);
    return 0;
}

int euler_study(int m, double t_final, const char *dossier)
{
    grid_t  g;
    int     n, *ia, *ja, err = 0, k;
    double *a, *evals, *evecs;
    double  alpha;
    char    nom_dat[256], nom_gp[256], png[256], commande[512];
    FILE   *f;

    MKDIR(dossier);

    if (grid_init(&g, m))
        return 1;

    if (prob_from_grid(&g, &n, &ia, &ja, &a)) {
        grid_free(&g);
        return 1;
    }

    evals = malloc(sizeof(double));
    evecs = malloc((size_t) n * sizeof(double));
    if (evals == NULL || evecs == NULL) {
        printf("\n ERREUR : pas assez de memoire pour le couple propre\n\n");
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    printf("\n=====================================================================\n");
    printf(" PROBLEME DE CAUCHY PAR EULER PROGRESSIF  --  m = %d  (n = %d)\n", m, n);
    printf("=====================================================================\n");

    if (primme(n, ia, ja, a, 1, evals, evecs)) {
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    alpha = critical_alpha(evals[0], REACTOR.beta2);
    printf("\n  Reference : lambda_min = %.10f m^-2, alpha = %.10f\n",
           evals[0], alpha);
    printf("  Dimension critique sur cette grille : L_c = %.6f m\n",
           alpha * REACTOR.L);
    printf("  Duree simulee : t_final = %g s\n", t_final);

    sprintf(nom_dat, "%s/euler_m%d.dat", dossier, m);
    sprintf(nom_gp,  "%s/euler_m%d.gp",  dossier, m);
    sprintf(png,     "%s/euler_m%d.png", dossier, m);

    f = fopen(nom_dat, "w");
    if (f == NULL) {
        printf("\n ERREUR : impossible d'ecrire '%s'\n\n", nom_dat);
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    for (k = 0; k < N_ECHELLES && !err; k++)
        err = integre(f, n, ia, ja, a, evecs, evals[0], alpha,
                      ECHELLES[k], t_final, k);

    fclose(f);

    if (!err) {
        err = write_script(nom_gp, nom_dat, png, m, t_final);
        if (!err) {
            sprintf(commande, "gnuplot \"%s\"", nom_gp);
            if (system(commande) != 0)
                printf("\n  ATTENTION : gnuplot n'a pas pu etre execute ;"
                       " lancez  gnuplot %s\n", nom_gp);
            else
                printf("\n  Figure ecrite : %s\n", png);
        }
    }

    free(ia); free(ja); free(a); free(evals); free(evecs);
    grid_free(&g);
    return err;
}

int euler_stability_demo(int m)
/*
   But
   ===
   Verifie experimentalement la borne de stabilite d'Euler progressif en
   integrant le cas critique avec des pas de temps encadrant

        dt_max = 2 / ( tau ( lambda_max - beta^2 ) ) .

   En dessous de cette borne le schema est stable ; au-dessus, le mode
   propre associe a lambda_max est amplifie a chaque pas et la solution
   explose, quelle que soit la physique du probleme.
*/
{
    grid_t  g;
    int     n, *ia, *ja, i, k, essai;
    double *a, *b, *evals, *evecs, *phi, *phi_np;
    double  alpha, facteur, lam_max, lam_gersh, dt_max, dt_gersh;
    double  norme, norme0, vmax;
    const double ratios[] = { 0.900, 0.990, 0.999, 1.001, 1.010, 1.100 };
    const int    n_essais = 6;

    if (grid_init(&g, m))
        return 1;

    if (prob_from_grid(&g, &n, &ia, &ja, &a)) {
        grid_free(&g);
        return 1;
    }

    evals = malloc(sizeof(double));
    evecs = malloc((size_t) n * sizeof(double));
    phi    = malloc((size_t) n * sizeof(double));
    phi_np = malloc((size_t) n * sizeof(double));

    if (!evals || !evecs || !phi || !phi_np
        || primme(n, ia, ja, a, 1, evals, evecs)) {
        free(ia); free(ja); free(a); free(evals); free(evecs);
        free(phi); free(phi_np); grid_free(&g);
        return 1;
    }

    /* mise a l'echelle critique : lambda_min devient exactement beta^2 */
    alpha   = critical_alpha(evals[0], REACTOR.beta2);
    facteur = 1.0 / (alpha * alpha);
    for (k = 0; k < ia[n]; k++)
        a[k] *= facteur;

    lam_gersh = gershgorin_max(n, ia, a);

    /* --- plus grande valeur propre exacte --------------------------------
       Le solveur ne fournit que les plus PETITES valeurs propres. On lui
       soumet donc B = sigma I - A avec sigma la borne de Gershgorin : B est
       symetrique, et sa plus petite valeur propre vaut sigma - lambda_max(A).
       Cette transformation ne coute aucune factorisation, seulement un
       changement de signe et un decalage de la diagonale. */
    b = malloc((size_t) ia[n] * sizeof(double));
    if (b == NULL) {
        printf("\n ERREUR : pas assez de memoire\n\n");
        free(ia); free(ja); free(a); free(evals); free(evecs);
        free(phi); free(phi_np); grid_free(&g);
        return 1;
    }
    for (i = 0; i < n; i++)
        for (k = ia[i]; k < ia[i + 1]; k++)
            b[k] = (ja[k] == i) ? lam_gersh - a[k] : -a[k];

    if (primme(n, ia, ja, b, 1, evals, evecs)) {
        free(b); free(ia); free(ja); free(a); free(evals); free(evecs);
        free(phi); free(phi_np); grid_free(&g);
        return 1;
    }
    lam_max = lam_gersh - evals[0];
    free(b);

    dt_gersh = 2.0 / (REACTOR.tau * (lam_gersh - REACTOR.beta2));
    dt_max   = 2.0 / (REACTOR.tau * (lam_max   - REACTOR.beta2));

    printf("\n=====================================================================\n");
    printf(" LIMITE DE STABILITE D'EULER PROGRESSIF  --  m = %d  (n = %d)\n", m, n);
    printf("=====================================================================\n");
    printf("\n  lambda_max exact                = %.6f m^-2\n", lam_max);
    printf("  lambda_max borne de Gershgorin  = %.6f m^-2   (majorant, +%.3f %%)\n",
           lam_gersh, 100.0 * (lam_gersh - lam_max) / lam_max);

    /* --- verification croisee par symetrie du spectre ---------------------
       Le graphe de la grille est biparti : en coloriant les noeuds selon la
       parite de i+j, aucune arete ne relie deux noeuds de meme couleur. Le
       spectre de la matrice d'adjacence est donc symetrique par rapport a 0.
       Comme A = (4/h^2) I - (1/h^2) Adj, le spectre de A est symetrique par
       rapport a 4/h^2, d'ou

            lambda_min + lambda_max  =  8/h^2  =  borne de Gershgorin .

       Cette identite relie deux calculs independants (la plus petite et la
       plus grande valeur propre) et les valide l'un par l'autre. */
    printf("\n  Verification par symetrie du spectre (graphe biparti) :\n");
    printf("    lambda_min + lambda_max = %.6f\n", REACTOR.beta2 + lam_max);
    printf("    8/h^2 (= borne de Gershgorin) = %.6f\n", lam_gersh);
    printf("    ecart = %.2e\n", fabs(REACTOR.beta2 + lam_max - lam_gersh));
    printf("\n  dt_max exact       = %.6e s\n", dt_max);
    printf("  dt_max Gershgorin  = %.6e s   (conservatif, c'est celui qu'on\n", dt_gersh);
    printf("                                    utilise en production : sur, et\n");
    printf("                                    obtenu sans calcul de valeur propre)\n\n");
    printf("  Les rapports ci-dessous sont relatifs au dt_max EXACT.\n\n");
    /* Le reacteur etant mis a l'echelle critique, la solution exacte tend
       vers le mode fondamental : || Phi || doit rester bornee. Toute
       croissance est donc imputable au schema, pas a la physique. C'est le
       critere utilise ici, bien plus sur qu'un simple test de depassement
       de capacite qui ne se declenche qu'apres des milliers de pas. */
    printf("   dt / dt_max   facteur |1-dt.tau(lmax-b^2)|   ||Phi||/||Phi(0)||   verdict\n");
    printf("  ---------------------------------------------------------------------------\n");

    norme0 = sqrt((double) n) * REACTOR.phi0;

    for (essai = 0; essai < n_essais; essai++) {
        double dt = ratios[essai] * dt_max;
        double g  = fabs(1.0 - dt * REACTOR.tau * (lam_max - REACTOR.beta2));

        for (i = 0; i < n; i++)
            phi[i] = REACTOR.phi0;

        norme = norme0;
        for (k = 0; k < N_PAS_DEMO; k++) {
            norme = euler_step(n, ia, ja, a, dt * REACTOR.tau, REACTOR.beta2,
                               phi, phi_np, &vmax);
            { double *tmp = phi; phi = phi_np; phi_np = tmp; }
            if (!(norme < 1e300)) break;
        }

        printf("      %5.3f              %10.6f            %12.4e     %s\n",
               ratios[essai], g, norme / norme0,
               (norme / norme0 < 10.0) ? "STABLE" : "DIVERGE");
    }

    printf("\n  -> la transition se produit exactement a dt = dt_max : le schema est\n");
    printf("     stable tant que le facteur d'amplification du mode le plus raide\n");
    printf("     reste <= 1, et diverge des qu'il le depasse.\n");
    printf("     La borne de Gershgorin majore lambda_max, donc le dt_max qu'elle\n");
    printf("     fournit est legerement trop petit : c'est exactement ce qu'on\n");
    printf("     attend d'un critere de securite, et il ne coute aucun calcul de\n");
    printf("     valeur propre. C'est celui utilise dans euler_study().\n\n");

    free(ia); free(ja); free(a); free(evals); free(evecs);
    free(phi); free(phi_np); grid_free(&g);
    return 0;
}
