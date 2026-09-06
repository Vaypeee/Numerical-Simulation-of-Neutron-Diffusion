#include <stdlib.h>
#include "primme.h"

/* variables statiques -- accessibles  */
static double *a;
static int n, *ia, *ja;

void matvec_primme(void *vx, void *vy, int *blockSize, primme_params *primme)
/*
   But
   ===
   Calcule le produit matrice-vecteur
                              vy = A*vx 
   pour le solveur aux valeurs propres PRIMME. La matrice A doit être
   "stoquée" au préalable dans les variables statiques 'n', 'ia', 'ja' et 'a'
   en utilisant le format CSR (Compressed Sparse Rows). Par "stoquer" 
   on veut dire ici stoquer la valeur de 'n' et les pointeurs vers les 
   tableaux 'ia', 'ja' et 'a'.

   Arguments
   =========
   vx        (input) - vecteur(s) d'entrée
   vy       (output) - vecteur(s) de produit A*vx
   blockSize (input) - nombre de vecteurs d'entrée
   primme    (input) - paramètres fournis par primme pour optimiser le calcul
                       (pas utilisé) 
*/
{
    int i, j, b;
    double *x = vx, *y=vy;

    for(b = 0; b < (*blockSize)*n; b+=n)
        for(i = 0; i < n; i++){
            y[b+i] = 0;
            for (j = ia[i]; j < ia[i + 1]; j++)
                y[b+i] += a[j] * x[b+ja[j]];
        }
} 

/**/

int primme(int primme_n, int *primme_ia, int *primme_ja, double *primme_a, 
           int nev, double *evals, double *evecs)

/*
   But
   ===
   Calcule les nev valeurs propres les plus basses de la matrice A de dimensions
   primme_n x primme_n stuquée dans le format CSR à l'aide des vecteurs
   primme_ia, primme_ja et primme_ia. 

  Arguments
  =========
  primme_n  (input) - le nombre d'inconnues dans le système
  primme_ia (input) - le tableau 'ia' de la matrice A
  primme_ja (input) - le tableau 'ja' de la matrice A
  primme_a  (input) - le tableau 'a' de la matrice A
  nev       (input) - le nombre de valeurs propres recherchées
  evals    (output) - le tableau des valeurs propres
  evecs    (output) - le tableau des vecteurs propres

  Retourne 0 si le calcul s'est bien déroulé, 1 si non.
*/
{
    int err;

    /* norme des résidus */
    double *resn = malloc(nev * sizeof(double));
    if (resn == NULL) {
        printf("\n ERREUR : pas assez de mémoire pour un vecteur auxilier dans la fonction primme\n\n");
        return 1;
    }

    /* sauvegarder les pointeurs dans des variables statiques */
    n = primme_n;
    a = primme_a;
    ja = primme_ja;
    ia = primme_ia;

    /* encoder les paramètres de PRIMME */
    primme_params primme;
    primme_initialize (&primme);
    primme.matrixMatvec = matvec_primme; 
                                /* nom de fonction du produit matrice-vecteur */
    primme.n = primme_n; /* dimensions de la matrice */
    primme.numEvals = nev; /* nombre de paires valeur propre-vecteur propre */
    primme.printLevel = 2; /* niveau d'affichage (1-4) */

    err = primme_set_method (DEFAULT_MIN_TIME, &primme);
    if(err){
        printf("\nPRIMME: erreur N %d dans le choix de la methode \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }
  
    /* afficher les papramètres de PRIMME */
    primme_display_params (primme);

    /* Caclul des valeurs et vecteurs propres */
    err = dprimme (evals, evecs, resn, &primme);
    if(err){
        printf("\nPRIMME: erreur N %d dans le calcul des valeurs propres \n    (voir 'Error Codes' dans le guide d'utilisateur)\n",err);
        return 1;
    }

    /* libérer la mémoire */
    primme_Free (&primme); free(resn);

    return 0;
}
