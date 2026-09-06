#ifndef RESIDUAL_H
#define RESIDUAL_H

/*
  Norme relative du residu d'un couple propre approche (lambda, x) :

                        || A x - lambda x ||_2
                        ----------------------
                              || x ||_2

  Deux implementations sont fournies : residual_norm(), fusionnee et sans
  allocation, et residual_norm_naive(), qui suit litteralement la formule.
  Elles doivent donner le meme resultat ; la seconde sert de reference de
  correction et de point de comparaison pour la mesure de performance.
*/

/* Implementation fusionnee : une seule passe sur la matrice, aucune allocation. */
double residual_norm(int n, const int *ia, const int *ja, const double *a,
                     double lambda, const double *x);

/* Implementation naive : produit matrice-vecteur dans un tableau temporaire,
   puis deux passes supplementaires pour les normes. Retourne -1 si
   l'allocation du tableau temporaire echoue. */
double residual_norm_naive(int n, const int *ia, const int *ja, const double *a,
                           double lambda, const double *x);

/* Compare les deux implementations en temps sur 'repet' evaluations et
   affiche le resultat. Sert a justifier le choix d'implementation. */
void residual_benchmark(int n, const int *ia, const int *ja, const double *a,
                        double lambda, const double *x, int repet);

#endif /* RESIDUAL_H */
