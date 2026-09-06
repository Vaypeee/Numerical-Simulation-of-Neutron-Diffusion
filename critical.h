#ifndef CRITICAL_H
#define CRITICAL_H

/*
  Dimensions critiques du reacteur (tache 3).

  Les dimensions fournies en annexe ne rendent pas le reacteur critique :
  la plus petite valeur propre de A n'est pas egale au buckling beta^2. On
  cherche donc le facteur d'homothetie alpha a appliquer a Omega pour que le
  reacteur soit stationnaire.
*/

/* Facteur d'homothetie rendant le reacteur critique : alpha = sqrt(lambda/beta2). */
double critical_alpha(double lambda_min, double beta2);

/* Affiche les dimensions du reacteur critique deduites de lambda_min. */
void   critical_report(double lambda_min, double h, int m);

/*
  Etude de convergence : resout le probleme aux valeurs propres pour une
  suite de grilles ou h est divise par deux a chaque etape, tabule
  lambda_min et la dimension critique, estime l'ordre de convergence
  observe et extrapole la limite h -> 0 par la methode de Richardson.
  Retourne 0 en cas de succes.
*/
int    critical_study(int niveaux);

#endif /* CRITICAL_H */
