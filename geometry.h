#ifndef GEOMETRY_H
#define GEOMETRY_H

/*
  Géométrie du réacteur et numérotation des inconnues.
  ====================================================

  Le domaine Ω est un carré [0,L] x [0,L] duquel on retire un rectangle
  aligné sur les axes :

        Ω  =  [0,L] x [0,L]  \  [cut_x0, cut_x1] x [cut_y0, cut_y1]

  Le rectangle retiré est FERMÉ : les nœuds situés sur son bord font donc
  partie de ∂Ω, où le flux s'annule ; ils ne sont pas des inconnues.

  Les inconnues sont les nœuds strictement intérieurs à Ω, numérotés dans
  l'ordre lexicographique : on part du coin sud-ouest, on parcourt la
  première rangée d'ouest en est, puis la deuxième, et ainsi de suite.
*/

/* Paramètres physiques et géométriques du réacteur (voir annexe de l'énoncé). */
typedef struct {
    double L;       /* côté du carré de départ                          [m]      */
    double cut_x0;  /* rectangle retiré : borne ouest                   [m]      */
    double cut_x1;  /* rectangle retiré : borne est                     [m]      */
    double cut_y0;  /* rectangle retiré : borne sud                     [m]      */
    double cut_y1;  /* rectangle retiré : borne nord                    [m]      */
    double beta2;   /* buckling matériel β²                             [m^-2]   */
    double tau;     /* paramètre cinétique τ                            [m² s^-1]*/
    double phi0;    /* flux initial constant φ₀                         [m^-2 s^-1] */
} reactor_t;

/* Le réacteur du projet, défini une seule fois dans geometry.c.
   Changer de numéro de projet = changer ces sept nombres, rien d'autre. */
extern const reactor_t REACTOR;

/* Grille de discrétisation et table de numérotation des inconnues. */
typedef struct {
    int     m;      /* nombre de points de grille par direction (bords compris) */
    double  h;      /* pas de discrétisation h = L/(m-1)                 [m]     */
    int     ci0;    /* indice i du bord ouest  du rectangle retiré               */
    int     ci1;    /* indice i du bord est    du rectangle retiré               */
    int     cj0;    /* indice j du bord sud    du rectangle retiré               */
    int     cj1;    /* indice j du bord nord   du rectangle retiré               */
    int     n;      /* nombre d'inconnues                                        */
    int    *num;    /* table m*m : numéro d'équation du nœud (i,j), -1 sinon     */
} grid_t;

/* Numéro d'équation du nœud (i,j) ; -1 si ce nœud n'est pas une inconnue. */
#define GRID_NUM(g, i, j)  ((g)->num[(i) + (g)->m * (j)])

/* Coordonnées physiques du nœud (i,j), en mètres. */
#define GRID_X(g, i)  ((i) * (g)->h)
#define GRID_Y(g, j)  ((j) * (g)->h)

/*
  Construit la grille à m points par direction et la table de numérotation.
  Retourne 0 en cas de succès, 1 en cas d'erreur (m invalide, grille
  incompatible avec la géométrie, ou mémoire insuffisante). En cas de
  grille incompatible, un message indique les valeurs de m admissibles.
*/
int  grid_init(grid_t *g, int m);

/* Libère la mémoire détenue par la grille. */
void grid_free(grid_t *g);

/* Le nœud (i,j) est-il une inconnue (c.-à-d. strictement intérieur à Ω) ? */
int  grid_is_unknown(const grid_t *g, int i, int j);

/* Le nœud (i,j) appartient-il à Ω, bord ∂Ω compris ? (utile pour l'affichage) */
int  grid_is_in_domain(const grid_t *g, int i, int j);

#endif /* GEOMETRY_H */
