#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "geometry.h"

/*
  Réacteur du projet n° 21 (voir annexe de l'énoncé) :

      Ω  =  [0, 7/2] x [0, 7/2]  \  [0, 2] x [0, 1]        (Figure 21)

  soit un carré de 3,5 m de côté amputé de son coin sud-ouest : une forme
  en L. Les paramètres physiques associés sont β² = 2 m^-2, τ = 100 m² s^-1
  et φ₀ = 4 m^-2 s^-1.
*/
const reactor_t REACTOR = {
    3.5,        /* L      [m]         */
    0.0, 2.0,   /* cut_x0, cut_x1 [m] */
    0.0, 1.0,   /* cut_y0, cut_y1 [m] */
    2.0,        /* β²     [m^-2]      */
    100.0,      /* τ      [m² s^-1]   */
    4.0         /* φ₀     [m^-2 s^-1] */
};

/* Tolérance relative utilisée pour décider qu'un quotient est entier. */
#define GEOM_TOL 1e-9

/*
  Convertit une abscisse (ou ordonnée) en indice de grille entier.
  Retourne -1 si la coordonnée ne tombe pas exactement sur un nœud, ce qui
  signifie que la grille ne respecte pas la géométrie du réacteur.
*/
static int coord_to_index(double v, int m)
{
    double t   = v * (m - 1) / REACTOR.L;
    int    idx = (int) floor(t + 0.5);

    return (fabs(t - idx) > GEOM_TOL * m) ? -1 : idx;
}

/* La grille à m points par direction est-elle compatible avec la géométrie ? */
static int grid_is_admissible(int m)
{
    return m >= 3
        && coord_to_index(REACTOR.cut_x0, m) >= 0
        && coord_to_index(REACTOR.cut_x1, m) >= 0
        && coord_to_index(REACTOR.cut_y0, m) >= 0
        && coord_to_index(REACTOR.cut_y1, m) >= 0;
}

/*
  Affiche les premières valeurs de m compatibles avec la géométrie.
  Sert uniquement à produire un message d'erreur utile.
*/
static void print_admissible_grids(int m_rejete)
{
    int m, trouves = 0;

    printf("\n ERREUR : m = %d est incompatible avec la geometrie du reacteur.\n", m_rejete);
    printf("          Les coins du rectangle retire doivent tomber sur des noeuds\n");
    printf("          de la grille, ce qui impose h = L/(m-1) diviseur exact de\n");
    printf("          %g m et %g m.\n", REACTOR.cut_x1, REACTOR.cut_y1);
    printf("          Valeurs de m admissibles :");

    for (m = 3; trouves < 8; m++) {
        if (grid_is_admissible(m)) {
            printf(" %d", m);
            trouves++;
        }
        if (m > 10000) break;   /* garde-fou : geometrie degeneree */
    }
    printf(" ...\n\n");
}

int grid_is_unknown(const grid_t *g, int i, int j)
{
    /* Hors du carré, ou sur son bord : condition de Dirichlet, pas une inconnue. */
    if (i <= 0 || i >= g->m - 1 || j <= 0 || j >= g->m - 1)
        return 0;

    /* Dans le rectangle retiré, bord compris : également sur ∂Ω. */
    if (i >= g->ci0 && i <= g->ci1 && j >= g->cj0 && j <= g->cj1)
        return 0;

    return 1;
}

int grid_is_in_domain(const grid_t *g, int i, int j)
{
    if (i < 0 || i > g->m - 1 || j < 0 || j > g->m - 1)
        return 0;

    /* Le bord du rectangle retiré appartient à Ω ; seul son intérieur en est exclu. */
    if (i > g->ci0 && i < g->ci1 && j > g->cj0 && j < g->cj1)
        return 0;

    return 1;
}

int grid_smallest_m(void)
{
    int m;

    for (m = 3; m <= 10000; m++)
        if (grid_is_admissible(m))
            return m;

    return 0;   /* geometrie degeneree : aucune grille ne convient */
}

int grid_init(grid_t *g, int m)
{
    int i, j, n = 0;

    if (m < 3) {
        printf("\n ERREUR : m = %d est trop petit (il faut m >= 3).\n\n", m);
        return 1;
    }

    if (!grid_is_admissible(m)) {
        print_admissible_grids(m);
        return 1;
    }

    g->m   = m;
    g->h   = REACTOR.L / (m - 1);
    g->ci0 = coord_to_index(REACTOR.cut_x0, m);
    g->ci1 = coord_to_index(REACTOR.cut_x1, m);
    g->cj0 = coord_to_index(REACTOR.cut_y0, m);
    g->cj1 = coord_to_index(REACTOR.cut_y1, m);

    g->num = malloc((size_t) m * m * sizeof(int));
    if (g->num == NULL) {
        printf("\n ERREUR : pas assez de memoire pour la table de numerotation\n\n");
        return 1;
    }

    /* Numérotation lexicographique : j (direction y) à l'extérieur, i à l'intérieur,
       en partant du coin sud-ouest. */
    for (j = 0; j < m; j++)
        for (i = 0; i < m; i++)
            GRID_NUM(g, i, j) = grid_is_unknown(g, i, j) ? n++ : -1;

    g->n = n;
    return 0;
}

void grid_free(grid_t *g)
{
    free(g->num);
    g->num = NULL;
    g->n   = 0;
}
