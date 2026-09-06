/* prototypes utilisés dans main.h */
double mytimer_cpu();
double mytimer_wall();
int prob(int m, int *n, int **ia, int **ja, double **a);
int primme(int primme_n, int* primme_ia, int* primme_ja, double* primme_a, 
           int nev, double *evals, double *evecs);

