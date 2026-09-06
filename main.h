#ifndef MAIN_H
#define MAIN_H

/*
  Pilote du projet. Les prototypes des différents modules se trouvent dans
  leurs en-têtes respectifs : geometry.h, prob.h, mytime.h, interface_primme.h.
*/

/* Numéro du projet, utilisé pour retrouver les fichiers CSR de référence
   ia.<N>.txt, ja.<N>.txt et a.<N>.txt fournis avec l'énoncé. */
#define PROJET_NUMERO 21

/* Grille de référence de l'annexe : celle dessinée en pointillé sur la figure. */
#define GRILLE_REFERENCE 8

/* Grilles utilisées pour la tâche 4 : l'énoncé demande un affichage du mode
   fondamental pour m <= 20 et pour m >= 500. Ces deux valeurs sont de la
   forme 7k+1, donc compatibles avec la géométrie du projet 21. */
#define PLOT_M_GROSSIER  15
#define PLOT_M_FIN      505

/* Tâche 5 : grille et durée par défaut de l'intégration en temps. La grille
   est un compromis entre finesse spatiale et coût, le pas de temps d'Euler
   progressif étant contraint en O(h²). */
#define EULER_M_DEFAUT   57
#define EULER_T_FINAL   1.0

#endif /* MAIN_H */
