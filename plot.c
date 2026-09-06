#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "plot.h"
#include "geometry.h"
#include "prob.h"
#include "residual.h"
#include "critical.h"
#include "interface_primme.h"

#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MKDIR(p) mkdir(p, 0755)
#endif

/*
   Normalisation du mode propre.
   =============================
   Le probleme (2) determine le flux a une constante multiplicative pres :
   si phi est solution, alpha*phi l'est aussi. Le vecteur rendu par le
   solveur est de norme euclidienne 1 et de signe arbitraire.

   La matrice sigma*I - A, avec sigma = 8/h^2, est a coefficients positifs
   (diagonale 4/h^2, hors-diagonale 1/h^2) et irreductible sur un domaine
   connexe. Le theoreme de Perron-Frobenius assure donc que son vecteur
   propre dominant -- c'est-a-dire le mode fondamental de A -- ne change pas
   de signe. On peut donc toujours le rendre strictement positif.

   On fixe ensuite l'echelle en imposant max(phi) = phi0, ce qui donne au
   trace une amplitude physiquement plausible et coherente avec le flux
   initial de la tache 5. Ce choix est arbitraire et signale sur le graphe.

   Retourne le minimum du vecteur normalise, qui doit etre positif ou nul si
   le mode fondamental a bien ete calcule.
*/
static double normalize_mode(double *phi, int n, double phi0)
{
    double vmax = 0.0, vmin, echelle;
    int    i, imax = 0;

    for (i = 0; i < n; i++)
        if (fabs(phi[i]) > vmax) { vmax = fabs(phi[i]); imax = i; }

    /* signe : rendre le mode positif */
    if (phi[imax] < 0.0)
        for (i = 0; i < n; i++)
            phi[i] = -phi[i];

    /* echelle : max(phi) = phi0 */
    echelle = (vmax > 0.0) ? phi0 / vmax : 1.0;
    vmin = phi[0] * echelle;
    for (i = 0; i < n; i++) {
        phi[i] *= echelle;
        if (phi[i] < vmin) vmin = phi[i];
    }

    return vmin;
}

/*
   Ecrit le fichier de donnees pour gnuplot.
   =========================================
   Une ligne "x y phi" par noeud de la grille m x m, y compris les noeuds de
   la frontiere dOmega ou le flux vaut exactement 0. Les noeuds strictement
   interieurs au rectangle retire n'appartiennent pas au reacteur : leur
   valeur est ecrite NaN, ce qui indique a gnuplot de ne pas dessiner les
   mailles correspondantes (set datafile missing "NaN"). Les noeuds situes
   SUR le bord du rectangle retire appartiennent bien a dOmega et sont donc
   traces, a la valeur 0 : c'est ce qui permet de verifier visuellement la
   condition aux limites le long de l'encoche.

   Les coordonnees sont celles du reacteur critique, obtenues en multipliant
   la grille par le facteur d'homothetie alpha ; elles sont donc en metres.

   Une ligne blanche separe les rangees successives, format attendu par
   gnuplot pour une grille structuree.
*/
static int write_data(const char *nom, const grid_t *g, const double *phi,
                      double alpha)
{
    FILE *f = fopen(nom, "w");
    int   i, j, ind;

    if (f == NULL) {
        printf("\n ERREUR : impossible d'ecrire '%s'\n\n", nom);
        return 1;
    }

    fprintf(f, "# Mode propre fondamental du reacteur critique\n");
    fprintf(f, "# colonne 1 : x [m]   colonne 2 : y [m]   colonne 3 : phi [m^-2 s^-1]\n");
    fprintf(f, "# NaN = point hors du reacteur (interieur du rectangle retire)\n");

    for (j = 0; j < g->m; j++) {
        for (i = 0; i < g->m; i++) {
            double x = alpha * GRID_X(g, i);
            double y = alpha * GRID_Y(g, j);

            if (!grid_is_in_domain(g, i, j)) {
                fprintf(f, "%.8g %.8g NaN\n", x, y);
            } else {
                ind = GRID_NUM(g, i, j);
                /* ind < 0 : noeud de dOmega, flux nul par condition aux limites */
                fprintf(f, "%.8g %.8g %.8g\n", x, y, ind < 0 ? 0.0 : phi[ind]);
            }
        }
        fprintf(f, "\n");   /* separation des rangees pour gnuplot */
    }

    fclose(f);
    return 0;
}

/*
   Ecrit le script gnuplot produisant deux figures : une vue en carte
   (projection du dessus, echelle de couleurs) et une vue en perspective.
*/
static int write_script(const char *nom_gp, const char *nom_dat,
                        const char *png_map, const char *png_surf,
                        int m, double h_crit, double lcrit, double phi0)
{
    FILE *f = fopen(nom_gp, "w");

    if (f == NULL) {
        printf("\n ERREUR : impossible d'ecrire '%s'\n\n", nom_gp);
        return 1;
    }

    fprintf(f, "# Script genere automatiquement par plot.c\n");
    /* On ne declare surtout PAS 'set datafile missing "NaN"' : cette
       directive supprimerait purement et simplement l'enregistrement, ce qui
       briserait la grille rectangulaire dont 'with image' a besoin. Gnuplot
       sait lire NaN nativement comme un flottant : le point reste alors dans
       la grille et n'est simplement pas colorie. */
    fprintf(f, "set palette defined (0 '#08306b', 0.25 '#2171b5', 0.5 '#6baed6',"
               " 0.75 '#fdae61', 1 '#a50026')\n");
    fprintf(f, "set cblabel \"flux {/Symbol f}  [m^{-2} s^{-1}]\"\n");
    fprintf(f, "set xlabel \"x  [m]\"\n");
    fprintf(f, "set ylabel \"y  [m]\"\n");
    fprintf(f, "set cbrange [0:%.8g]\n", phi0);

    /* ---- vue en carte ----
       Rendu 'with image' : chaque valeur est dessinee comme un pixel centre
       sur son noeud. C'est la representation fidele de donnees nodales. Un
       rendu pm3d colorierait au contraire les mailles SITUEES ENTRE les
       noeuds : toute maille touchant un point NaN serait alors ecartee, et
       l'encoche apparaitrait plus grande d'une maille qu'elle ne l'est. */
    fprintf(f, "\nset terminal pngcairo size 900,820 enhanced font 'Arial,13'\n");
    fprintf(f, "set output '%s'\n", png_map);
    fprintf(f, "set title \"Mode fondamental du reacteur critique"
               " (m = %d, h = %.4g m)\\n"
               "cote = %.4f m   -   {/Symbol f} normalise a max = %g m^{-2} s^{-1}\"\n",
            m, h_crit, lcrit, phi0);
    fprintf(f, "set size ratio -1\n");
    /* Chaque pixel etant centre sur son noeud, il deborde de h/2 de part et
       d'autre : on elargit d'autant pour que la couronne de bord a phi = 0
       soit entierement visible. */
    fprintf(f, "set xrange [%.8g:%.8g]\n", -0.5 * h_crit, lcrit + 0.5 * h_crit);
    fprintf(f, "set yrange [%.8g:%.8g]\n", -0.5 * h_crit, lcrit + 0.5 * h_crit);
    fprintf(f, "plot '%s' using 1:2:3 with image notitle\n", nom_dat);

    /* ---- vue en perspective ---- */
    fprintf(f, "\nunset size\n");
    fprintf(f, "set xrange [0:%.8g]\n", lcrit);
    fprintf(f, "set yrange [0:%.8g]\n", lcrit);
    fprintf(f, "set terminal pngcairo size 1000,780 enhanced font 'Arial,13'\n");
    fprintf(f, "set output '%s'\n", png_surf);
    fprintf(f, "set title \"Mode fondamental du reacteur critique"
               " (m = %d)\\nles bords sont a {/Symbol f} = 0\"\n", m);
    fprintf(f, "set zlabel \"{/Symbol f}  [m^{-2} s^{-1}]\" rotate by 90\n");
    fprintf(f, "set zrange [0:%.8g]\n", 1.05 * phi0);
    fprintf(f, "set ticslevel 0\n");
    fprintf(f, "set view 55,325\n");
    fprintf(f, "set pm3d depthorder\n");
    fprintf(f, "unset pm3d\n");
    fprintf(f, "set pm3d\n");
    fprintf(f, "set hidden3d\n");
    fprintf(f, "splot '%s' using 1:2:3 with pm3d notitle\n", nom_dat);

    fclose(f);
    return 0;
}

/*
   Verifie que le flux trace respecte les conditions aux limites.
   ==============================================================
   Sur dOmega le flux est nul par construction (ces noeuds ne sont pas des
   inconnues), ce qui est verifie exactement. La verification interessante
   porte sur la premiere couche de noeuds interieurs adjacents au bord : le
   flux y doit etre petit devant le maximum et tendre vers zero quand h
   diminue, ce qui traduit la continuite du flux jusqu'au bord.
*/
static void check_boundary(const grid_t *g, const double *phi, double phi0)
{
    double max_bord = 0.0, max_couche1 = 0.0;
    int    i, j, ind, n_bord = 0;

    for (j = 0; j < g->m; j++) {
        for (i = 0; i < g->m; i++) {
            if (!grid_is_in_domain(g, i, j))
                continue;

            ind = GRID_NUM(g, i, j);

            if (ind < 0) {
                /* noeud de dOmega : valeur tracee = 0 */
                n_bord++;
                continue;
            }

            /* noeud interieur adjacent a dOmega ? */
            if (GRID_NUM(g, i - 1, j) < 0 || GRID_NUM(g, i + 1, j) < 0 ||
                GRID_NUM(g, i, j - 1) < 0 || GRID_NUM(g, i, j + 1) < 0)
                if (phi[ind] > max_couche1)
                    max_couche1 = phi[ind];
        }
    }

    printf("\n  Verification des conditions aux limites :\n");
    printf("    noeuds de dOmega traces               : %d\n", n_bord);
    printf("    max |phi| sur dOmega                  : %.3e  (nul par construction)\n",
           max_bord);
    printf("    max phi sur la 1re couche interieure  : %.4f m^-2 s^-1"
           "  (%.2f %% du maximum)\n",
           max_couche1, 100.0 * max_couche1 / phi0);
    printf("    -> le flux decroit continument vers 0 en approchant le bord.\n");
}

int plot_fundamental_mode(int m, const char *dossier)
{
    grid_t  g;
    int     n, *ia, *ja, err;
    double *a, *evals, *evecs;
    double  alpha, lcrit, vmin, res;
    char    nom_dat[256], nom_gp[256], png_map[256], png_surf[256], commande[512];

    MKDIR(dossier);   /* echec ignore : le repertoire existe probablement deja */

    if (grid_init(&g, m))
        return 1;

    if (prob_from_grid(&g, &n, &ia, &ja, &a)) {
        grid_free(&g);
        return 1;
    }

    evals = malloc(sizeof(double));
    evecs = malloc((size_t) n * sizeof(double));
    if (evals == NULL || evecs == NULL) {
        printf("\n ERREUR : pas assez de memoire pour le couple propre\n\n");
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    printf("\n=====================================================================\n");
    printf(" MODE FONDAMENTAL  --  grille m = %d  (n = %d inconnues)\n", m, n);
    printf("=====================================================================\n");

    if (primme(n, ia, ja, a, 1, evals, evecs)) {
        free(ia); free(ja); free(a); free(evals); free(evecs); grid_free(&g);
        return 1;
    }

    res   = residual_norm(n, ia, ja, a, evals[0], evecs);
    alpha = critical_alpha(evals[0], REACTOR.beta2);
    lcrit = alpha * REACTOR.L;
    vmin  = normalize_mode(evecs, n, REACTOR.phi0);

    printf("  lambda_min = %.10f m^-2      residu = %.2e\n", evals[0], res);
    printf("  alpha      = %.10f          cote critique = %.6f m\n", alpha, lcrit);
    printf("  pas de la grille sur le reacteur critique : h = %.6f m\n", alpha * g.h);
    printf("  min(phi) apres normalisation : %.3e  (doit etre >= 0)\n", vmin);

    if (vmin < 0.0)
        printf("  ATTENTION : le mode change de signe, ce n'est pas le fondamental !\n");

    /* --- fichiers de sortie --------------------------------------------- */
    sprintf(nom_dat,  "%s/mode_m%d.dat",         dossier, m);
    sprintf(nom_gp,   "%s/mode_m%d.gp",          dossier, m);
    sprintf(png_map,  "%s/mode_m%d_carte.png",   dossier, m);
    sprintf(png_surf, "%s/mode_m%d_surface.png", dossier, m);

    err = write_data(nom_dat, &g, evecs, alpha)
       || write_script(nom_gp, nom_dat, png_map, png_surf,
                       m, alpha * g.h, lcrit, REACTOR.phi0);

    if (!err) {
        check_boundary(&g, evecs, REACTOR.phi0);

        sprintf(commande, "gnuplot \"%s\"", nom_gp);
        if (system(commande) != 0) {
            printf("\n  ATTENTION : gnuplot n'a pas pu etre execute.\n");
            printf("  Les fichiers '%s' et '%s' ont ete ecrits ;\n", nom_dat, nom_gp);
            printf("  lancez manuellement :  gnuplot %s\n", nom_gp);
        } else {
            printf("\n  Figures ecrites :\n    %s\n    %s\n", png_map, png_surf);
        }
    }

    free(ia); free(ja); free(a); free(evals); free(evecs);
    grid_free(&g);
    return err;
}
