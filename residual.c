#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "residual.h"
#include "mytime.h"

double residual_norm(int n, const int *ia, const int *ja, const double *a,
                     double lambda, const double *x)
/*
   But
   ===
   Evalue la norme relative du residu

                || A x - lambda x ||_2  /  || x ||_2

   pour un couple (lambda, x) approchant une paire propre de A, la matrice
   etant stockee au format CSR par 'n', 'ia', 'ja' et 'a'.

   Strategie d'implementation
   ==========================
   Le calcul naif enchaine trois etapes : y = A x, puis || y - lambda x ||,
   puis || x ||. Il alloue un vecteur temporaire y de taille n, l'ecrit une
   fois puis le relit, et relit x deux fois de plus. Le cout en trafic
   memoire est donc de nnz + 5n lectures/ecritures de doubles en plus de la
   matrice elle-meme.

   Ici tout est fusionne en une seule boucle sur les lignes :

     - la ligne i de A n'est parcourue qu'une fois, et sa contribution est
       immediatement combinee avec -lambda*x[i] ;
     - le residu de la ligne, r, ne quitte jamais les registres : il n'est
       jamais ecrit en memoire, et aucun tableau temporaire n'est alloue ;
     - x[i] est charge une seule fois dans la variable xi, qui sert a la
       fois pour le terme -lambda*x[i] et pour la norme de x ;
     - les bornes de ligne ia[i] et ia[i+1] sont conservees d'une iteration
       a l'autre (k_fin devient k_deb), ce qui evite une relecture de ia ;
     - les deux sommes de carres sont accumulees dans des variables locales,
       donc dans des registres.

   Le trafic memoire se reduit ainsi a une lecture unique de 'ja', de 'a',
   de 'ia' et de 'x' : c'est le minimum possible, puisque chacune de ces
   donnees est indispensable au resultat.

   Le mot-cle 'restrict' indique au compilateur que les quatre tableaux ne
   se recouvrent pas, ce qui lui permet de garder les accumulateurs en
   registre et de vectoriser la boucle interne.

   Arguments
   =========
   n      (input) - nombre d'inconnues
   ia     (input) - tableau 'ia' de la matrice A
   ja     (input) - tableau 'ja' de la matrice A
   a      (input) - tableau 'a' de la matrice A
   lambda (input) - valeur propre approchee
   x      (input) - vecteur propre approche

   Sortie
   ======
   La norme relative du residu, ou 0 si x est le vecteur nul.
*/
{
    const int    * const restrict pia = ia;
    const int    * const restrict pja = ja;
    const double * const restrict pa  = a;
    const double * const restrict px  = x;

    double norme_r2 = 0.0;   /* || A x - lambda x ||^2 */
    double norme_x2 = 0.0;   /* || x ||^2              */
    int    i, k;
    int    k_deb = pia[0];

    for (i = 0; i < n; i++) {
        const int    k_fin = pia[i + 1];
        const double xi    = px[i];
        double       r     = -lambda * xi;   /* reste dans un registre */

        for (k = k_deb; k < k_fin; k++)
            r += pa[k] * px[pja[k]];

        norme_r2 += r * r;
        norme_x2 += xi * xi;

        k_deb = k_fin;                       /* evite de relire ia[i+1] */
    }

    return (norme_x2 > 0.0) ? sqrt(norme_r2 / norme_x2) : 0.0;
}

double residual_norm_naive(int n, const int *ia, const int *ja, const double *a,
                           double lambda, const double *x)
/*
   But
   ===
   Meme calcul que residual_norm, mais en suivant litteralement la formule :
   un produit matrice-vecteur stocke dans un tableau temporaire, puis les
   deux normes. Sert de reference de correction et de temoin pour la mesure
   de performance.

   Sortie
   ======
   La norme relative du residu, ou -1 si l'allocation echoue.
*/
{
    double *y = malloc((size_t) n * sizeof(double));
    double  norme_r2 = 0.0, norme_x2 = 0.0, d;
    int     i, k;

    if (y == NULL)
        return -1.0;

    /* etape 1 : y = A x */
    for (i = 0; i < n; i++) {
        y[i] = 0.0;
        for (k = ia[i]; k < ia[i + 1]; k++)
            y[i] += a[k] * x[ja[k]];
    }

    /* etape 2 : || y - lambda x ||^2 */
    for (i = 0; i < n; i++) {
        d = y[i] - lambda * x[i];
        norme_r2 += d * d;
    }

    /* etape 3 : || x ||^2 */
    for (i = 0; i < n; i++)
        norme_x2 += x[i] * x[i];

    free(y);
    return (norme_x2 > 0.0) ? sqrt(norme_r2 / norme_x2) : 0.0;
}

void residual_benchmark(int n, const int *ia, const int *ja, const double *a,
                        double lambda, const double *x, int repet)
/*
   Mesure le temps des deux implementations sur 'repet' evaluations et
   affiche la comparaison. La valeur retournee par chaque appel est
   accumulee dans une variable puis affichee, afin que le compilateur ne
   puisse pas eliminer les appels comme du code mort.
*/
{
    double t1, t2, t_fus, t_naif, s_fus = 0.0, s_naif = 0.0;
    int    k;

    t1 = mytimer_wall();
    for (k = 0; k < repet; k++)
        s_fus += residual_norm(n, ia, ja, a, lambda, x);
    t2 = mytimer_wall();
    t_fus = t2 - t1;

    t1 = mytimer_wall();
    for (k = 0; k < repet; k++)
        s_naif += residual_norm_naive(n, ia, ja, a, lambda, x);
    t2 = mytimer_wall();
    t_naif = t2 - t1;

    printf("\nCOMPARAISON DES IMPLEMENTATIONS DU RESIDU (%d evaluations) :\n", repet);
    printf("  fusionnee (1 passe, 0 allocation) : %8.4f s\n", t_fus);
    printf("  naive     (3 passes + tableau y)  : %8.4f s\n", t_naif);
    if (t_fus > 0.0)
        printf("  acceleration                      : %8.2fx\n", t_naif / t_fus);
    printf("  ecart sur les valeurs calculees   : %.3e\n",
           fabs(s_fus - s_naif) / (repet > 0 ? repet : 1));
}
