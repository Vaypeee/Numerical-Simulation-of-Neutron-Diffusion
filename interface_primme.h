#ifndef INTERFACE_PRIMME_H
#define INTERFACE_PRIMME_H

/* Regle la verbosite du solveur : 0 = silencieux (defaut), 2 = trace des
   iterations, 3 = trace + parametres complets. */
void primme_set_verbosity(int v);

/* Compteur de produits matrice-vecteur effectues par le solveur. */
void primme_reset_matvec_count(void);
long primme_get_matvec_count(void);

/* Regle la tolerance de convergence (critere ||r|| < eps ||A||).
   eps = 0 retablit le defaut de PRIMME. */
void primme_set_tolerance(double eps);

/* prototype */
int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs);

#endif /* INTERFACE_PRIMME_H */
