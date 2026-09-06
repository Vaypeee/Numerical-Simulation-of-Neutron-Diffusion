#ifndef PROB_H
#define PROB_H

#include "geometry.h"

/*
  Génère la matrice A du problème aux valeurs propres A Φ = β² Φ issu de la
  discrétisation par différences finies centrées de l'opérateur de Laplace
  sur la géométrie du réacteur, au format CSR.

  Voir prob.c pour la description détaillée des arguments.
*/
int prob(int m, int *n, int **ia, int **ja, double **a);

/*
  Variante qui réutilise une grille déjà construite et la laisse à
  l'appelant (utile pour l'affichage, qui a besoin de la numérotation).
*/
int prob_from_grid(const grid_t *g, int *n, int **ia, int **ja, double **a);

#endif /* PROB_H */
