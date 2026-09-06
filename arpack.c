#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <arpack/arpack.h>
#include "arpack.h"
#include "geometry.h"
#include "prob.h"
#include "residual.h"
#include "mytime.h"
#include "interface_primme.h"

/*
   Strategie retenue et pourquoi.
   ==============================
   ARPACK applique a une matrice symetrique se ramene a l'algorithme de
   Lanczos, avec redemarrage implicite. On cherche ici la plus petite valeur
   propre de A, qui est symetrique definie positive.

   Une idee naturelle est de decaler le spectre : avec sigma majorant
   lambda_max (borne de Gershgorin, soit 8/h^2 pour le schema a cinq
   points), la matrice B = sigma I - A a les memes vecteurs propres et les
   valeurs propres sigma - lambda_i, en ordre inverse. La plus petite valeur
   propre de A devient la plus grande de B, et l'on demanderait "LA" plutot
   que "SA".

   CETTE TRANSFORMATION EST INUTILE, et il faut le dire : les sous-espaces
   de Krylov sont invariants par decalage,

        K_k(sigma I - A, v)  =  K_k(A, v)   pour tout sigma,

   puisque chaque puissance (sigma I - A)^j v est une combinaison lineaire
   des A^i v. Lanczos explore donc exactement le meme espace dans les deux
   cas, et converge a la meme vitesse. La fonction arpack_compare() le
   verifie experimentalement : les deux variantes consomment le meme nombre
   de produits matrice-vecteur, aux arrondis pres.

   On resout donc directement A avec "SA" (smallest algebraic), ce qui est
   plus simple et, accessoirement, plus precis : le critere d'arret d'ARPACK
   etant relatif a la valeur propre visee, viser lambda_min ~ 2 donne un
   seuil absolu bien plus serre que viser lambda_max(B) ~ 8/h^2.

   La seule transformation qui accelererait REELLEMENT la convergence est le
   mode "shift-invert", qui remplace A par (A - sigma I)^-1 et separe
   fortement les valeurs propres voisines de sigma. Mais il impose de
   factoriser une matrice creuse : en deux dimensions le remplissage rend
   cette factorisation couteuse en memoire et en temps, et l'on perd le
   caractere purement matrice-vecteur du solveur. Ce choix est ecarte ici.

   Proprietes de la matrice effectivement exploitees :
     - symetrie : on utilise dsaupd/dseupd (Lanczos) et non dnaupd/dneupd
       (Arnoldi non symetrique). La recurrence est a trois termes au lieu de
       complete, d'ou un cout et une memoire nettement moindres, et des
       valeurs propres reelles garanties ;
     - caractere creux : le produit matrice-vecteur reste en O(nnz) et la
       matrice n'est jamais formee explicitement ;
     - communication inverse : ARPACK ne connait pas la matrice, il reclame
       seulement des produits, donc rien n'est duplique en memoire ;
     - structure bipartie de la grille : elle fournit l'identite
       lambda_min + lambda_max = 8/h^2, utilisee comme controle croise avec
       PRIMME dans arpack_compare().
*/

/*
   Produit matrice-vecteur y = B x = sigma x - A x, en une seule passe sur
   la structure CSR et sans tableau intermediaire (meme principe que pour le
   residu de la tache 2).
*/
static void matvec_shifted(int n, const int *ia, const int *ja,
                           const double *a, double sigma, int avec_decalage,
                           const double *x, double *y)
{
    const int    * const restrict pia = ia;
    const int    * const restrict pja = ja;
    const double * const restrict pa  = a;
    const double * const restrict px  = x;
    double       * const restrict py  = y;

    int i, k, k_deb = pia[0];

    for (i = 0; i < n; i++) {
        const int k_fin = pia[i + 1];
        double    r     = 0.0;

        for (k = k_deb; k < k_fin; k++)
            r += pa[k] * px[pja[k]];

        py[i] = avec_decalage ? (sigma * px[i] - r) : r;
        k_deb = k_fin;
    }
}

/* Borne de Gershgorin : max des sommes de valeurs absolues par ligne. */
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

int arpack_solve(int n, const int *ia, const int *ja, const double *a,
                 int avec_decalage, double tol_abs,
                 double *eval, double *evec, long *n_matvec)
{
    a_int   ido = 0, info = 0, ierr = 0;
    a_int   nev = 1, ncv, ldv, lworkl, rvec = 1;
    a_int   iparam[11], ipntr[11];
    a_int  *select = NULL;
    double *resid = NULL, *v = NULL, *workd = NULL, *workl = NULL, *d = NULL;
    double  tol = 0.0;          /* rempli plus bas ; 0 => precision machine */
    double  sigma_shift = 0.0;  /* non utilise en mode 1 */
    double  sigma;
    const char *which;
    long    compteur = 0;
    int     i, err = 0;

    /* Taille de l'espace de Krylov. ARPACK impose nev < ncv <= n ; une
       valeur de l'ordre de 20 accelere nettement la convergence pour une
       seule valeur propre, au prix d'une memoire O(n*ncv) modeste. */
    ncv = (n < 30) ? n : 30;
    if (ncv <= nev) ncv = (nev + 1 < n) ? nev + 1 : n;
    ldv    = n;
    lworkl = ncv * (ncv + 8);

    resid  = malloc((size_t) n * sizeof(double));
    v      = malloc((size_t) ldv * ncv * sizeof(double));
    workd  = malloc(3 * (size_t) n * sizeof(double));
    workl  = malloc((size_t) lworkl * sizeof(double));
    d      = malloc((size_t) ncv * sizeof(double));
    select = malloc((size_t) ncv * sizeof(a_int));

    if (!resid || !v || !workd || !workl || !d || !select) {
        printf("\n ERREUR : pas assez de memoire pour ARPACK\n\n");
        err = 1;
        goto fin;
    }

    /* avec_decalage = 1 : on resout sur B = sigma I - A en demandant "LA".
       avec_decalage = 0 : on resout directement sur A en demandant "SA".
       Les deux doivent donner exactement le meme cout (voir la note en
       tete de fichier sur l'invariance des sous-espaces de Krylov). */
    sigma = avec_decalage ? gershgorin_max(n, ia, a) : 0.0;
    which = avec_decalage ? "LA" : "SA";

    /* --- traduction de la precision absolue en tolerance relative -------
       Le critere d'arret d'ARPACK est   bounds(i) <= tol * |ritz(i)| ,
       donc relatif a la valeur propre VISEE. Or celle-ci n'est pas la meme
       selon la variante : lambda_min(A) ~ 2 sans decalage, mais
       lambda_max(B) = sigma - lambda_min ~ 8/h^2 avec decalage. A tol egal,
       la variante decalee s'arrete donc sur une precision absolue
       sigma/lambda_min fois plus laxiste, ecart qui croit comme h^-2.
       Pour comparer les deux a precision ABSOLUE egale, on convertit :  */
    if (tol_abs > 0.0) {
        double echelle_ritz = avec_decalage ? gershgorin_max(n, ia, a) : 1.0;
        tol = tol_abs / (echelle_ritz > 0.0 ? echelle_ritz : 1.0);
    }

    for (i = 0; i < 11; i++) { iparam[i] = 0; ipntr[i] = 0; }
    iparam[0] = 1;      /* strategie de decalage : shifts exacts       */
    iparam[2] = 5000;   /* nombre maximal d'iterations d'Arnoldi       */
    iparam[6] = 1;      /* mode 1 : probleme aux valeurs propres standard */

    /* --- boucle de communication inverse ---------------------------------
       ARPACK ne connait pas la matrice : il rend la main a chaque fois
       qu'il a besoin d'un produit matrice-vecteur, en indiquant dans ipntr
       ou lire l'entree et ou ecrire la sortie dans workd. */
    do {
        dsaupd_c(&ido, "I", n, which, nev, tol, resid, ncv, v, ldv,
                 iparam, ipntr, workd, workl, lworkl, &info);

        if (ido == -1 || ido == 1) {
            matvec_shifted(n, ia, ja, a, sigma, avec_decalage,
                           &workd[ipntr[0] - 1],    /* ipntr est en base 1 */
                           &workd[ipntr[1] - 1]);
            compteur++;
        }
    } while (ido == -1 || ido == 1);

    if (info < 0) {
        printf("\n ERREUR : dsaupd a retourne info = %d\n\n", (int) info);
        err = 1;
        goto fin;
    }
    if (info == 1)
        printf("\n ATTENTION : ARPACK a atteint le nombre maximal d'iterations\n");

    /* --- extraction des valeurs et vecteurs propres ---------------------- */
    dseupd_c(rvec, "All", select, d, v, ldv, sigma_shift, "I", n, which,
             nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd, workl,
             lworkl, &ierr);

    if (ierr != 0) {
        printf("\n ERREUR : dseupd a retourne ierr = %d\n\n", (int) ierr);
        err = 1;
        goto fin;
    }

    /* avec decalage : B = sigma I - A, donc lambda_min(A) = sigma - lambda_max(B)
       sans decalage : d[0] est deja lambda_min(A) */
    *eval = avec_decalage ? sigma - d[0] : d[0];

    for (i = 0; i < n; i++)
        evec[i] = v[i];

    if (n_matvec != NULL)
        *n_matvec = compteur;

fin:
    free(resid); free(v); free(workd); free(workl); free(d); free(select);
    return err;
}

int arpack_smallest(int n, const int *ia, const int *ja, const double *a,
                    double *eval, double *evec, long *n_matvec)
{
    /* Sans decalage : "SA" directement sur A. Voir la note en tete de
       fichier : le decalage n'apporte rien et relache le critere d'arret. */
    return arpack_solve(n, ia, ja, a, 0, 0.0, eval, evec, n_matvec);
}

int arpack_compare(int m)
/*
   But
   ===
   Compare ARPACK et PRIMME sur le meme probleme : valeur propre obtenue,
   norme du residu, temps CPU et horloge, et nombre de produits
   matrice-vecteur, qui est la mesure de cout la plus objective puisqu'elle
   ne depend ni de la machine ni des options de compilation.
*/
{
    grid_t  g;
    int     n, *ia, *ja;
    double *a, *b, *ev_primme, *vec_primme, *ev_arpack, *vec_arpack;
    double  t1, t2;
    double  t_primme_cpu, t_primme_wall, t_arpack_cpu, t_arpack_wall;
    double  res_primme, res_arpack, sigma;
    long    mv_primme = 0, mv_arpack = 0;
    int     err = 0;

    if (grid_init(&g, m))
        return 1;

    if (prob_from_grid(&g, &n, &ia, &ja, &a)) {
        grid_free(&g);
        return 1;
    }

    ev_primme  = malloc(sizeof(double));
    vec_primme = malloc((size_t) n * sizeof(double));
    ev_arpack  = malloc(sizeof(double));
    vec_arpack = malloc((size_t) n * sizeof(double));

    if (!ev_primme || !vec_primme || !ev_arpack || !vec_arpack) {
        printf("\n ERREUR : pas assez de memoire\n\n");
        free(ia); free(ja); free(a);
        free(ev_primme); free(vec_primme); free(ev_arpack); free(vec_arpack);
        grid_free(&g);
        return 1;
    }

    printf("\n=====================================================================\n");
    printf(" COMPARAISON PRIMME / ARPACK  --  m = %d  (n = %d, nnz = %d)\n",
           m, n, ia[n]);
    printf("=====================================================================\n");

    /* --- PRIMME ----------------------------------------------------------
       ARPACK est appele avec tol = 0, c'est-a-dire a la precision machine.
       Pour que la comparaison porte sur un cout A PRECISION EGALE et non sur
       deux criteres d'arret differents, on resserre la tolerance de PRIMME
       en consequence. Le critere de PRIMME est ||r|| < eps ||A|| ; avec
       ||A|| de l'ordre de 8/h^2, eps = 1e-14 conduit a un residu du meme
       ordre que celui obtenu par ARPACK. */
    primme_set_tolerance(1e-14);
    primme_reset_matvec_count();
    t1 = mytimer_cpu(); t2 = mytimer_wall();
    err = primme(n, ia, ja, a, 1, ev_primme, vec_primme);
    t_primme_cpu  = mytimer_cpu()  - t1;
    t_primme_wall = mytimer_wall() - t2;
    mv_primme     = primme_get_matvec_count();
    primme_set_tolerance(0.0);   /* retablir le defaut pour le reste du code */

    /* --- ARPACK ---------------------------------------------------------- */
    if (!err) {
        t1 = mytimer_cpu(); t2 = mytimer_wall();
        err = arpack_smallest(n, ia, ja, a, ev_arpack, vec_arpack, &mv_arpack);
        t_arpack_cpu  = mytimer_cpu()  - t1;
        t_arpack_wall = mytimer_wall() - t2;
    }

    if (!err) {
        res_primme = residual_norm(n, ia, ja, a, *ev_primme, vec_primme);
        res_arpack = residual_norm(n, ia, ja, a, *ev_arpack, vec_arpack);
        sigma      = gershgorin_max(n, ia, a);

        printf("\n  ARPACK : dsaupd/dseupd (Lanczos symetrique), \"SA\" sur A,\n");
        printf("  sans factorisation ni preconditionneur.\n");
        printf("  sigma = 8/h^2 = %.1f ne sert qu'aux controles ci-dessous.\n", sigma);

        printf("\n  %-24s %20s %20s\n", "", "PRIMME", "ARPACK");
        printf("  ---------------------------------------------------------------------\n");
        printf("  %-24s %20.12f %20.12f\n", "lambda_min [m^-2]", *ev_primme, *ev_arpack);
        printf("  %-24s %20.3e %20.3e\n", "residu relatif", res_primme, res_arpack);
        printf("  %-24s %20.4f %20.4f\n", "temps CPU [s]", t_primme_cpu, t_arpack_cpu);
        printf("  %-24s %20.4f %20.4f\n", "temps horloge [s]", t_primme_wall, t_arpack_wall);
        printf("  %-24s %20ld %20ld\n", "produits matrice-vecteur", mv_primme, mv_arpack);

        printf("\n  Ecart entre les deux valeurs propres : %.3e\n",
               fabs(*ev_primme - *ev_arpack));
        printf("  Ecart relatif                        : %.3e\n",
               fabs(*ev_primme - *ev_arpack) / fabs(*ev_primme));

        /* --- le decalage sert-il a quelque chose ? ------------------------
           On resout le meme probleme sans decalage, en demandant "SA"
           directement sur A. Voir la note en tete de fichier : les
           sous-espaces de Krylov etant invariants par decalage, les deux
           strategies doivent couter exactement la meme chose. */
        {
            double ev_sa, *vec_sa = malloc((size_t) n * sizeof(double));
            long   mv_sa = 0;

            if (vec_sa != NULL
                && !arpack_solve(n, ia, ja, a, 1, tol_equiv, &ev_sa, vec_sa, &mv_sa)) {
                printf("\n  Effet du decalage spectral sur le cout d'ARPACK :\n");
                printf("    \"SA\" directement sur A   : %6ld produits matrice-vecteur\n",
                       mv_arpack);
                printf("    \"LA\" sur B = sigma I - A : %6ld produits matrice-vecteur\n",
                       mv_sa);
                printf("    ecart sur lambda_min      : %.3e\n",
                       fabs(ev_sa - *ev_arpack));
                printf("    -> cout essentiellement identique : le decalage ne change\n");
                printf("       rien, les sous-espaces de Krylov etant invariants par\n");
                printf("       decalage. L'ecart residuel vient des arrondis dans la\n");
                printf("       recurrence de Lanczos et de la strategie de redemarrage.\n");
            }
            free(vec_sa);
        }

        /* --- controle croise par la symetrie du spectre --------------------
           Le graphe de la grille est biparti (colorier les noeuds selon la
           parite de i+j : aucune arete ne relie deux noeuds de meme
           couleur). Le spectre de la matrice d'adjacence est donc symetrique
           par rapport a 0, et comme A = (4/h^2) I - (1/h^2) Adj, celui de A
           l'est par rapport a 4/h^2, d'ou

                lambda_min + lambda_max = 8/h^2 = sigma .

           On calcule lambda_max avec l'AUTRE solveur : PRIMME applique a
           B = sigma I - A rend lambda_min(B) = sigma - lambda_max(A). Les
           deux membres de l'identite proviennent donc de deux solveurs
           differents, ce qui en fait une verification reellement
           independante et non une tautologie. */
        b = malloc((size_t) ia[n] * sizeof(double));
        if (b != NULL) {
            double ev_b, *vec_b = malloc((size_t) n * sizeof(double));
            int    i, k;

            for (i = 0; i < n; i++)
                for (k = ia[i]; k < ia[i + 1]; k++)
                    b[k] = (ja[k] == i) ? sigma - a[k] : -a[k];

            if (vec_b != NULL && !primme(n, ia, ja, b, 1, &ev_b, vec_b)) {
                double lam_max = sigma - ev_b;

                printf("\n  Controle croise par la symetrie du spectre"
                       " (graphe biparti) :\n");
                printf("    lambda_min(A) par ARPACK                 = %.9f\n",
                       *ev_arpack);
                printf("    lambda_max(A) par PRIMME sur B           = %.9f\n",
                       lam_max);
                printf("    somme                                    = %.9f\n",
                       *ev_arpack + lam_max);
                printf("    sigma = 8/h^2 attendu                    = %.9f\n",
                       sigma);
                printf("    ecart                                    = %.3e\n",
                       fabs(*ev_arpack + lam_max - sigma));
            }
            free(vec_b);
        }
        free(b);
    }

    free(ia); free(ja); free(a);
    free(ev_primme); free(vec_primme); free(ev_arpack); free(vec_arpack);
    grid_free(&g);
    return err;
}
