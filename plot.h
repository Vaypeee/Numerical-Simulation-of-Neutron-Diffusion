#ifndef PLOT_H
#define PLOT_H

/*
  Visualisation du mode propre fondamental (tache 4).

  Le flux est trace sur le reacteur CRITIQUE, c'est-a-dire sur le domaine
  dilate alpha*Omega : c'est celui-la qui est stationnaire. Les points de la
  frontiere dOmega sont inclus dans le trace, avec un flux nul.

  Le mode propre n'etant defini qu'a une constante multiplicative pres, il
  est normalise de facon a etre positif et a valoir phi0 en son maximum.
*/

/*
  Calcule le mode fondamental sur la grille a m points, ecrit les fichiers
  de donnees et de commandes gnuplot, lance gnuplot et verifie les
  conditions aux limites.

  m      (input) - nombre de points de grille par direction
  dossier (input) - repertoire de sortie (cree s'il n'existe pas)

  Retourne 0 en cas de succes, 1 en cas d'erreur.
*/
int plot_fundamental_mode(int m, const char *dossier);

#endif /* PLOT_H */
