#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "main.h"
#include "geometry.h"
#include "prob.h"
#include "mytime.h"
#include "interface_primme.h"

/*
  Lit tous les nombres réels séparés par des blancs contenus dans un fichier.
  Retourne le tableau alloué (à libérer par l'appelant) et écrit le nombre de
  valeurs lues dans *count ; retourne NULL si le fichier est illisible.
*/
static double *read_numbers(const char *filename, int *count)
{
    FILE   *f = fopen(filename, "r");
    double *v = NULL, *tmp;
    int     cap = 0, k = 0;
    double  x;

    if (f == NULL) {
        printf(" ERREUR : impossible d'ouvrir '%s'\n", filename);
        return NULL;
    }

    while (fscanf(f, "%lf", &x) == 1) {
        if (k == cap) {
            cap = cap ? 2 * cap : 256;
            tmp = realloc(v, (size_t) cap * sizeof(double));
            if (tmp == NULL) { free(v); fclose(f); return NULL; }
            v = tmp;
        }
        v[k++] = x;
    }

    fclose(f);
    *count = k;
    return v;
}

/*
  Compare la matrice CSR produite par prob() aux fichiers de référence
  ia.<N>.txt, ja.<N>.txt et a.<N>.txt fournis avec l'énoncé.
  Retourne 0 si tout concorde, 1 sinon.
*/
static int check_reference(int n, const int *ia, const int *ja, const double *a)
{
    char    f_ia[64], f_ja[64], f_a[64];
    double *r_ia = NULL, *r_ja = NULL, *r_a = NULL;
    int     n_ia = 0, n_ja = 0, n_a = 0;
    int     k, erreurs = 0, nnz = ia[n];

    sprintf(f_ia, "ia.%d.txt", PROJET_NUMERO);
    sprintf(f_ja, "ja.%d.txt", PROJET_NUMERO);
    sprintf(f_a,  "a.%d.txt",  PROJET_NUMERO);

    r_ia = read_numbers(f_ia, &n_ia);
    r_ja = read_numbers(f_ja, &n_ja);
    r_a  = read_numbers(f_a,  &n_a);

    if (r_ia == NULL || r_ja == NULL || r_a == NULL) {
        free(r_ia); free(r_ja); free(r_a);
        return 1;
    }

    printf("\nVERIFICATION par rapport aux fichiers CSR de reference :\n");
    printf("  reference : n = %d, nnz = %d\n", n_ia - 1, n_ja);
    printf("  calcule   : n = %d, nnz = %d\n", n, nnz);

    if (n_ia - 1 != n || n_ja != nnz || n_a != nnz) {
        printf("  -> DIMENSIONS DIFFERENTES\n\n");
        free(r_ia); free(r_ja); free(r_a);
        return 1;
    }

    for (k = 0; k <= n; k++)
        if (ia[k] != (int) r_ia[k]) {
            if (erreurs < 5)
                printf("  ia[%d] : calcule %d, reference %d\n", k, ia[k], (int) r_ia[k]);
            erreurs++;
        }

    for (k = 0; k < nnz; k++) {
        if (ja[k] != (int) r_ja[k]) {
            if (erreurs < 5)
                printf("  ja[%d] : calcule %d, reference %d\n", k, ja[k], (int) r_ja[k]);
            erreurs++;
        }
        if (fabs(a[k] - r_a[k]) > 1e-12 * fabs(r_a[k])) {
            if (erreurs < 5)
                printf("  a[%d]  : calcule %.15g, reference %.15g\n", k, a[k], r_a[k]);
            erreurs++;
        }
    }

    if (erreurs == 0)
        printf("  -> IDENTIQUE : ia, ja et a concordent sur %d elements non nuls.\n\n", nnz);
    else
        printf("  -> %d DIFFERENCE(S) detectee(s).\n\n", erreurs);

    free(r_ia); free(r_ja); free(r_a);
    return erreurs != 0;
}

static void usage(const char *prog)
{
    printf("\nUsage : %s [-m M] [--check] [--no-solve]\n", prog);
    printf("  -m M        nombre de points de grille par direction (defaut : %d)\n",
           GRILLE_REFERENCE);
    printf("  --check     comparer la matrice aux fichiers CSR de reference\n");
    printf("  --no-solve  ne pas appeler PRIMME (generation de la matrice seule)\n\n");
}

int main(int argc, char *argv[])
{
    int     m = GRILLE_REFERENCE, nev = 1, k;
    int     do_check = 0, do_solve = 1;
    int     n, *ia, *ja;
    double *a, *evals, *evecs;
    double  tc1, tc2, tw1, tw2;
    grid_t  g;

    /* --- lecture des arguments ------------------------------------------ */
    for (k = 1; k < argc; k++) {
        if (strcmp(argv[k], "-m") == 0 && k + 1 < argc)      m = atoi(argv[++k]);
        else if (strcmp(argv[k], "--check") == 0)            do_check = 1;
        else if (strcmp(argv[k], "--no-solve") == 0)         do_solve = 0;
        else { usage(argv[0]); return 1; }
    }

    /* --- generation du probleme ----------------------------------------- */
    if (grid_init(&g, m))
        return 1;

    if (prob_from_grid(&g, &n, &ia, &ja, &a)) {
        grid_free(&g);
        return 1;
    }

    printf("\nREACTEUR : Omega = [0,%g]x[0,%g] \\ [%g,%g]x[%g,%g]  (projet n %d)\n",
           REACTOR.L, REACTOR.L, REACTOR.cut_x0, REACTOR.cut_x1,
           REACTOR.cut_y0, REACTOR.cut_y1, PROJET_NUMERO);
    printf("           beta^2 = %g m^-2   tau = %g m^2/s   phi0 = %g m^-2 s^-1\n",
           REACTOR.beta2, REACTOR.tau, REACTOR.phi0);
    printf("\nPROBLEME : m = %5d   h = %10.6f m   n = %8d   nnz = %9d\n",
           m, g.h, n, ia[n]);

    if (do_check && check_reference(n, ia, ja, a)) {
        free(ia); free(ja); free(a); grid_free(&g);
        return 1;
    }

    if (!do_solve) {
        free(ia); free(ja); free(a); grid_free(&g);
        return 0;
    }

    /* --- allouer la memoire pour vecteurs & valeurs propres -------------- */
    evals = malloc((size_t) nev * sizeof(double));
    evecs = malloc((size_t) nev * n * sizeof(double));

    if (evals == NULL || evecs == NULL) {
        printf("\n ERREUR : pas assez de memoire pour les vecteurs et valeurs propres\n\n");
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    /* --- primme : resolution -------------------------------------------- */
    tc1 = mytimer_cpu(); tw1 = mytimer_wall();
    if (primme(n, ia, ja, a, nev, evals, evecs)) {
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }
    tc2 = mytimer_cpu(); tw2 = mytimer_wall();

    printf("\nTemps de solution (CPU)     : %8.3f s", tc2 - tc1);
    printf("\nTemps de solution (horloge) : %8.3f s\n", tw2 - tw1);
    printf("\nValeur propre minimale calculee : %.12g m^-2\n\n", evals[0]);

    free(ia); free(ja); free(a); free(evals); free(evecs);
    grid_free(&g);
    return 0;
}
