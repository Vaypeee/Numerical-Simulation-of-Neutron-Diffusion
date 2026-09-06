#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "critical.h"
#include "geometry.h"
#include "prob.h"
#include "residual.h"
#include "interface_primme.h"

/*
   Pourquoi une homothetie de rapport alpha divise-t-elle A par alpha^2 ?
   =====================================================================
   Soit Omega' = alpha*Omega le domaine dilate. En posant x' = alpha*x, une
   fonction phi definie sur Omega se transporte en phi'(x') = phi(x'/alpha).
   Chaque derivation par rapport a x' fait apparaitre un facteur 1/alpha par
   la regle de derivation en chaine, donc le laplacien se transforme selon

        Laplacien_{x'} phi'  =  (1/alpha^2) (Laplacien_x phi) .

   Au niveau discret, la meme chose se lit directement sur la matrice : si on
   garde le meme nombre m de points par direction, le pas devient
   h' = alpha*h, et comme tous les coefficients de A sont proportionnels a
   1/h^2 (4/h^2 sur la diagonale, -1/h^2 hors diagonale),

        A(alpha*Omega)  =  A(Omega) / alpha^2 .

   La structure creuse et la numerotation sont inchangees : seule l'echelle
   des coefficients varie. Les valeurs propres suivent donc

        lambda(alpha)  =  lambda(1) / alpha^2 ,

   avec les memes vecteurs propres. La condition de criticite
   lambda_min(alpha) = beta^2 donne alors immediatement

        alpha  =  sqrt( lambda_min(1) / beta^2 ) .
*/

double critical_alpha(double lambda_min, double beta2)
{
    return sqrt(lambda_min / beta2);
}

/*
   Nombre de coins rentrants du domaine.
   =====================================
   Un coin du rectangle retire situe strictement a l'interieur du carre cree
   un coin rentrant : le domaine y presente un angle interieur de 3*pi/2. Un
   coin situe sur le bord du carre, lui, ne fait que decouper la frontiere
   sans creer de singularite.
*/
static int count_reentrant_corners(void)
{
    const double cx[2] = { REACTOR.cut_x0, REACTOR.cut_x1 };
    const double cy[2] = { REACTOR.cut_y0, REACTOR.cut_y1 };
    int i, j, compte = 0;

    for (i = 0; i < 2; i++)
        for (j = 0; j < 2; j++)
            if (cx[i] > 0.0 && cx[i] < REACTOR.L && cy[j] > 0.0 && cy[j] < REACTOR.L)
                compte++;

    return compte;
}

/*
   Ordre de convergence theorique attendu sur la valeur propre.
   ============================================================
   Sur un domaine convexe, le schema a cinq points est d'ordre 2 et la
   valeur propre converge en O(h^2).

   Des qu'il existe un coin rentrant d'angle interieur omega = 3*pi/2, la
   fonction propre n'est plus reguliere : au voisinage du coin elle se
   comporte comme r^(pi/omega) = r^(2/3), ou r est la distance au coin. La
   derivee seconde y est singuliere, l'erreur de troncature du schema n'est
   plus O(h^2), et l'erreur sur la valeur propre se degrade en

        O( h^(2*pi/omega) )  =  O( h^(4/3) ) .

   Pour la geometrie du projet 21, le coin (2,1) est le seul coin rentrant,
   et on attend donc p = 4/3 ~= 1.333.
*/
static double theoretical_order(void)
{
    return (count_reentrant_corners() > 0) ? 4.0 / 3.0 : 2.0;
}

void critical_report(double lambda_min, double h, int m)
/*
   Affiche les dimensions du reacteur critique deduites de lambda_min,
   calcule sur la grille a m points de pas h.
*/
{
    double alpha = critical_alpha(lambda_min, REACTOR.beta2);
    double aire  = (REACTOR.L * REACTOR.L
                    - (REACTOR.cut_x1 - REACTOR.cut_x0)
                    * (REACTOR.cut_y1 - REACTOR.cut_y0));

    printf("\nDIMENSIONS CRITIQUES (grille m = %d, h = %g m) :\n", m, h);
    printf("  lambda_min              = %.10f m^-2\n", lambda_min);
    printf("  beta^2                  = %.10f m^-2\n", REACTOR.beta2);
    printf("  alpha = sqrt(lambda_min/beta^2) = %.10f\n", alpha);
    printf("  le reacteur donne est %s : %s\n",
           lambda_min > REACTOR.beta2 ? "SOUS-CRITIQUE" : "SUR-CRITIQUE",
           lambda_min > REACTOR.beta2 ? "il faut l'agrandir." : "il faut le reduire.");
    printf("\n  Geometrie critique  Omega_c = alpha * Omega :\n");
    printf("    cote du carre         : %.7f m   (au lieu de %g m)\n",
           alpha * REACTOR.L, REACTOR.L);
    printf("    rectangle retire      : [%.7f, %.7f] x [%.7f, %.7f] m\n",
           alpha * REACTOR.cut_x0, alpha * REACTOR.cut_x1,
           alpha * REACTOR.cut_y0, alpha * REACTOR.cut_y1);
    printf("    aire du coeur         : %.7f m^2\n", alpha * alpha * aire);
}

/*
   Resout le probleme aux valeurs propres sur la grille a m points et
   retourne la plus petite valeur propre dans *lambda, le pas dans *h, le
   nombre d'inconnues dans *n_out et le residu relatif dans *res_out.
   Retourne 0 en cas de succes.
*/
static int solve_min_eigenvalue(int m, double *lambda, double *h, int *n_out,
                                double *res_out)
{
    grid_t  g;
    int     n, *ia, *ja, err;
    double *a, *evals, *evecs;

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

    err = primme(n, ia, ja, a, 1, evals, evecs);

    if (!err) {
        *lambda  = evals[0];
        *h       = g.h;
        *n_out   = n;
        *res_out = residual_norm(n, ia, ja, a, evals[0], evecs);
    }

    free(ia); free(ja); free(a); free(evals); free(evecs);
    grid_free(&g);
    return err;
}

int critical_study(int niveaux)
/*
   But
   ===
   Tabule la dimension critique du reacteur en fonction du pas de
   discretisation, sur une suite de grilles ou h est divise par deux a
   chaque etape : m_k = (m_0 - 1) * 2^k + 1, ou m_0 est la grille de
   reference de l'annexe. Comme h est divise par deux exactement, l'ordre de
   convergence observe se lit directement :

       si  q(h) = q* + C h^p ,   alors   p = log2( (q_k - q_k+1)
                                                 / (q_k+1 - q_k+2) )

   et l'extrapolation de Richardson donne l'estimation amelioree

       q*  ~=  q_n  +  (q_n - q_n-1) / (2^p - 1) .

   L'ordre est estime sur la dimension critique L_c = alpha*L, qui est la
   grandeur physique demandee.

   Arguments
   =========
   niveaux (input) - nombre de grilles de la suite

   Sortie
   ======
   0 - execution avec succes
   1 - erreurs
*/
{
    int    *m_tab, *n_tab;
    double *lam, *lc, *h_tab, *res_tab;
    int     m0 = grid_smallest_m();
    int     k, err = 0, p_compte = 0;
    double  p_moy = 0.0, richardson;

    if (m0 == 0) {
        printf("\n ERREUR : aucune grille admissible pour cette geometrie\n\n");
        return 1;
    }

    m_tab   = malloc((size_t) niveaux * sizeof(int));
    n_tab   = malloc((size_t) niveaux * sizeof(int));
    lam     = malloc((size_t) niveaux * sizeof(double));
    lc      = malloc((size_t) niveaux * sizeof(double));
    h_tab   = malloc((size_t) niveaux * sizeof(double));
    res_tab = malloc((size_t) niveaux * sizeof(double));

    if (!m_tab || !n_tab || !lam || !lc || !h_tab || !res_tab) {
        printf("\n ERREUR : pas assez de memoire pour l'etude de convergence\n\n");
        free(m_tab); free(n_tab); free(lam); free(lc); free(h_tab); free(res_tab);
        return 1;
    }

    printf("\n=====================================================================\n");
    printf(" CONVERGENCE DE LA DIMENSION CRITIQUE EN FONCTION DE h\n");
    printf(" (raffinement dyadique : h divise par 2 a chaque ligne)\n");
    printf("=====================================================================\n\n");
    printf("    m        h [m]        n    lambda_min [m^-2]   L_crit [m]   residu\n");
    printf("  --------------------------------------------------------------------\n");

    for (k = 0; k < niveaux; k++) {
        m_tab[k] = (m0 - 1) * (1 << k) + 1;

        if (solve_min_eigenvalue(m_tab[k], &lam[k], &h_tab[k], &n_tab[k], &res_tab[k])) {
            err = 1;
            break;
        }

        lc[k] = critical_alpha(lam[k], REACTOR.beta2) * REACTOR.L;

        printf("  %5d  %10.7f %8d  %16.10f  %11.7f  %7.1e\n",
               m_tab[k], h_tab[k], n_tab[k], lam[k], lc[k], res_tab[k]);
        fflush(stdout);
    }

    if (!err && niveaux >= 3) {
        printf("\n  Variations successives de L_crit et ordre observe :\n\n");
        printf("    m=%-5d -> m=%-5d : delta = %11.3e\n",
               m_tab[0], m_tab[1], lc[1] - lc[0]);

        for (k = 1; k < niveaux - 1; k++) {
            double d1 = lc[k]     - lc[k - 1];
            double d2 = lc[k + 1] - lc[k];

            printf("    m=%-5d -> m=%-5d : delta = %11.3e",
                   m_tab[k], m_tab[k + 1], d2);

            /* L'ordre n'a de sens que si deux variations consecutives sont de
               meme signe ; au debut de la suite le regime asymptotique n'est
               pas encore atteint et la suite change de sens. */
            if (d1 * d2 > 0.0 && d2 != 0.0) {
                /* Le dernier ordre calcule est le plus proche du regime
                   asymptotique : c'est celui qu'on retient pour Richardson. */
                p_moy = log(fabs(d1 / d2)) / log(2.0);
                p_compte++;
                printf("   ordre p = %6.3f", p_moy);
            } else {
                printf("   (regime asymptotique non atteint)");
            }
            printf("\n");
        }

        if (p_compte > 0) {
            double p_theo = theoretical_order();
            double rich_theo;

            k          = niveaux - 1;
            richardson = lc[k] + (lc[k] - lc[k - 1]) / (pow(2.0, p_moy) - 1.0);
            rich_theo  = lc[k] + (lc[k] - lc[k - 1]) / (pow(2.0, p_theo) - 1.0);

            printf("\n  Coins rentrants du domaine : %d\n", count_reentrant_corners());
            printf("  Ordre theorique attendu    : p = %.3f", p_theo);
            printf(count_reentrant_corners() > 0
                   ? "   (coin rentrant : phi ~ r^(2/3))\n"
                   : "   (domaine convexe : schema d'ordre 2)\n");
            printf("  Ordre observe (dernier)    : p = %.3f\n", p_moy);

            printf("\n  Extrapolation de Richardson sur L_crit :\n");
            printf("    avec p observe  = %.3f : L_crit(h -> 0) = %.7f m\n", p_moy, richardson);
            printf("    avec p theorique = %.3f : L_crit(h -> 0) = %.7f m\n", p_theo, rich_theo);
            printf("    ecart entre les deux estimations       : %.2e m\n",
                   fabs(richardson - rich_theo));
            printf("  Erreur restante sur la grille la plus fine : %.2e m  (%.4f %%)\n",
                   fabs(richardson - lc[k]),
                   100.0 * fabs(richardson - lc[k]) / richardson);

            printf("\n  REACTEUR CRITIQUE (valeur retenue : %.4f m) :\n", richardson);
            printf("    cote du carre    : %.6f m\n", richardson);
            printf("    rectangle retire : [%.6f, %.6f] x [%.6f, %.6f] m\n",
                   richardson / REACTOR.L * REACTOR.cut_x0,
                   richardson / REACTOR.L * REACTOR.cut_x1,
                   richardson / REACTOR.L * REACTOR.cut_y0,
                   richardson / REACTOR.L * REACTOR.cut_y1);
        }
        printf("\n");
    }

    free(m_tab); free(n_tab); free(lam); free(lc); free(h_tab); free(res_tab);
    return err;
}
