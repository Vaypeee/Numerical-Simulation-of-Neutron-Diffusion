#ifndef ARPACK_H
#define ARPACK_H

/*
  Solveur aux valeurs propres alternatif : ARPACK (tache 6).

  ARPACK (Arnoldi Package, licence BSD, librement disponible sur
  https://github.com/opencollab/arpack-ng) implemente la methode d'Arnoldi
  redemarree implicitement. Pour une matrice symetrique elle se reduit a
  l'algorithme de Lanczos, qui n'utilise qu'une recurrence a trois termes.

  Le solveur est utilise ici en exploitant au maximum les proprietes de la
  matrice : symetrie, caractere creux, definie positivite, et structure
  bipartie de la grille. Voir arpack.c pour le detail.
*/

/*
  Calcule la plus petite valeur propre de A et le vecteur propre associe.

  n, ia, ja, a (input)  - matrice au format CSR
  eval        (output)  - plus petite valeur propre
  evec        (output)  - vecteur propre associe (tableau de taille n)
  n_matvec    (output)  - nombre de produits matrice-vecteur effectues
                          (peut etre NULL)

  Retourne 0 en cas de succes, 1 en cas d'erreur.
*/
int arpack_smallest(int n, const int *ia, const int *ja, const double *a,
                    double *eval, double *evec, long *n_matvec);

/*
  Variante permettant de choisir la strategie :
    avec_decalage = 1 : resout sur B = sigma I - A en demandant "LA"
    avec_decalage = 0 : resout directement sur A en demandant "SA"
  tol_rel est la tolerance RELATIVE transmise a ARPACK (critere
  bounds <= tol * |ritz|) ; 0 demande la precision machine. Comme la
  valeur propre visee differe d'une variante a l'autre, c'est a
  l'appelant de convertir s'il veut une precision absolue identique.
  Sert a comparer les deux strategies a precision egale.
*/
int arpack_solve(int n, const int *ia, const int *ja, const double *a,
                 int avec_decalage, double tol_rel,
                 double *eval, double *evec, long *n_matvec);

/*
  Etude comparative complete PRIMME / ARPACK sur la grille a m points :
  valeur propre, residu, temps de calcul et nombre de produits
  matrice-vecteur. Retourne 0 en cas de succes.
*/
int arpack_compare(int m);

#endif /* ARPACK_H */
