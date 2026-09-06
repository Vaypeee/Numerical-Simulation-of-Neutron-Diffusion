#include <stdlib.h>
#include <stdio.h>
#include "prob.h"
#include "geometry.h"

int prob_from_grid(const grid_t *g, int *n, int **ia, int **ja, double **a)
/*
   But
   ===
   Génère la matrice n x n qui correspond à la discrétisation sur une grille
   cartésienne régulière de l'opérateur de Laplace à deux dimensions

            d    d        d    d
         - == ( == u ) - == ( == u )        sur Omega
           dx   dx       dy   dy

   avec la fonction u qui satisfait les conditions aux limites de Dirichlet
   u = 0 sur dOmega. Le domaine Omega est le carré [0,L] x [0,L] amputé d'un
   rectangle, tel que décrit dans geometry.h.

   La discrétisation utilise la formule de différence centrée pour la dérivée
   seconde, ce qui donne en chaque nœud intérieur (i,j) :

       -(phi[i-1][j] - 2 phi[i][j] + phi[i+1][j]) / h^2
       -(phi[i][j-1] - 2 phi[i][j] + phi[i][j+1]) / h^2  =  beta^2 phi[i][j]

   Un voisin situé sur dOmega a un flux nul : le terme correspondant
   disparaît de l'équation, mais la contribution -2 phi[i][j] du nœud
   central subsiste. L'élément diagonal vaut donc 4/h^2 en tout nœud, quel
   que soit le nombre de voisins éliminés.

   La numérotation des inconnues est lexicographique, la direction x étant
   parcourue avant celle de y ; seuls les nœuds intérieurs à Omega sont des
   inconnues. Comme la numérotation est croissante dans cet ordre, les
   voisins sud et ouest ont toujours un numéro inférieur à celui du nœud
   central, et les voisins est et nord un numéro supérieur : en remplissant
   la ligne dans l'ordre sud, ouest, diagonale, est, nord, les indices de
   colonne sont automatiquement triés par ordre croissant.

   La matrice est retournée dans le format CSR qui est défini par le
   scalaire 'n' et les trois tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   g  (input)  - grille de discrétisation, déjà construite par grid_init
   n  (output) - pointeur vers le nombre d'inconnues dans le système
   ia (output) - pointeur vers le tableau 'ia' de la matrice A
   ja (output) - pointeur vers le tableau 'ja' de la matrice A
   a  (output) - pointeur vers le tableau 'a' de la matrice A

   Sortie
   ======
   0 - exécution avec succès
   1 - erreurs
*/
{
    int     i, j, ind, nnz = 0;
    int    *ja_tmp;
    double *a_tmp;
    double  invh2 = 1.0 / (g->h * g->h);

    *n = g->n;

    /* Chaque ligne compte au plus 5 éléments non nuls : on alloue cette borne
       supérieure, puis on rétrécit les tableaux une fois nnz connu. */
    *ia = malloc(((size_t) g->n + 1) * sizeof(int));
    *ja = malloc(5 * (size_t) g->n * sizeof(int));
    *a  = malloc(5 * (size_t) g->n * sizeof(double));

    if (*ia == NULL || *ja == NULL || *a == NULL) {
        printf("\n ERREUR : pas assez de memoire pour generer la matrice\n\n");
        free(*ia); free(*ja); free(*a);
        *ia = NULL; *ja = NULL; *a = NULL;
        return 1;
    }

    /* Partie principale : remplissage de la matrice, ligne par ligne, dans
       l'ordre lexicographique des nœuds intérieurs. */
    for (j = 0; j < g->m; j++) {
        for (i = 0; i < g->m; i++) {

            ind = GRID_NUM(g, i, j);
            if (ind < 0)
                continue;           /* nœud sur dOmega : pas d'équation */

            /* marquer le début de la ligne dans le tableau 'ia' */
            (*ia)[ind] = nnz;

            /* voisin sud */
            if (GRID_NUM(g, i, j - 1) >= 0) {
                (*a)[nnz]  = -invh2;
                (*ja)[nnz] = GRID_NUM(g, i, j - 1);
                nnz++;
            }

            /* voisin ouest */
            if (GRID_NUM(g, i - 1, j) >= 0) {
                (*a)[nnz]  = -invh2;
                (*ja)[nnz] = GRID_NUM(g, i - 1, j);
                nnz++;
            }

            /* élément diagonal */
            (*a)[nnz]  = 4.0 * invh2;
            (*ja)[nnz] = ind;
            nnz++;

            /* voisin est */
            if (GRID_NUM(g, i + 1, j) >= 0) {
                (*a)[nnz]  = -invh2;
                (*ja)[nnz] = GRID_NUM(g, i + 1, j);
                nnz++;
            }

            /* voisin nord */
            if (GRID_NUM(g, i, j + 1) >= 0) {
                (*a)[nnz]  = -invh2;
                (*ja)[nnz] = GRID_NUM(g, i, j + 1);
                nnz++;
            }
        }
    }

    /* dernier élément du tableau 'ia' */
    (*ia)[g->n] = nnz;

    /* rétrécir 'ja' et 'a' à la taille effectivement utilisée */
    ja_tmp = realloc(*ja, (size_t) nnz * sizeof(int));
    a_tmp  = realloc(*a,  (size_t) nnz * sizeof(double));
    if (ja_tmp != NULL) *ja = ja_tmp;
    if (a_tmp  != NULL) *a  = a_tmp;

    return 0;
}

int prob(int m, int *n, int **ia, int **ja, double **a)
/*
   But
   ===
   Idem que prob_from_grid, mais construit et libère la grille en interne.
   C'est la signature attendue par main.c et l'interface de PRIMME.

   Arguments
   =========
   m (input)   - nombre de points par direction dans la grille
   n  (output) - pointeur vers le nombre d'inconnues dans le système
   ia (output) - pointeur vers le tableau 'ia' de la matrice A
   ja (output) - pointeur vers le tableau 'ja' de la matrice A
   a  (output) - pointeur vers le tableau 'a' de la matrice A

   Sortie
   ======
   0 - exécution avec succès
   1 - erreurs
*/
{
    grid_t g;
    int    err;

    if (grid_init(&g, m))
        return 1;

    err = prob_from_grid(&g, n, ia, ja, a);
    grid_free(&g);

    return err;
}
