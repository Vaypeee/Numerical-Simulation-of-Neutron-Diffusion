#ifndef EULER_H
#define EULER_H

/*
  Resolution du probleme de Cauchy vectoriel par la methode d'Euler
  progressive (tache 5).

      dPhi/dt  =  -tau ( A Phi(t) - beta^2 Phi(t) ) ,   t dans [0, +inf[
      Phi(0)   =  1 * phi0

  On compare l'evolution pour des reacteurs legerement plus petits, de la
  taille critique exacte, et legerement plus grands que la dimension
  critique determinee a la tache 3.
*/

/*
  Etude complete : integre le probleme de Cauchy pour plusieurs facteurs
  d'echelle autour de la dimension critique, ecrit les donnees et produit la
  figure gnuplot de l'evolution de || Phi(t) ||.

  m       (input) - nombre de points de grille par direction
  t_final (input) - duree simulee [s]
  dossier (input) - repertoire de sortie

  Retourne 0 en cas de succes.
*/
int euler_study(int m, double t_final, const char *dossier);

/*
  Demonstration de la limite de stabilite : integre le meme probleme avec
  des pas de temps situes de part et d'autre de la borne
  dt_max = 2 / (tau (lambda_max - beta^2)) et montre que le schema diverge
  des qu'elle est franchie.

  Retourne 0 en cas de succes.
*/
int euler_stability_demo(int m);

#endif /* EULER_H */
