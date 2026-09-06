# =====================================================================
#  Projet de complements de calcul numerique -- projet n 21
#  Equation de diffusion neutronique 2D stationnaire
# =====================================================================

# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/

# BLAS et LAPACK. OpenBLAS fournit les deux ; sur une distribution qui
# livre les deux bibliotheques separement, remplacer par -lblas -llapack.
LIBBLAS = -lopenblas

# ARPACK : solveur aux valeurs propres alternatif (tache 6), licence BSD
LIBARPACK = -larpack

# toutes les librairies
LIB = $(LIBP) $(LIBARPACK) $(LIBBLAS) -lm

COPT = -O3 -Wall -Wextra

OBJ = main.o geometry.o prob.o residual.o critical.o plot.o euler.o arpack.o mytime.o interface_primme.o

default: main

main: $(OBJ)
	$(CC) $(COPT) $^ -o $@ $(LIB)

main.o: main.c main.h geometry.h prob.h residual.h critical.h plot.h euler.h arpack.h mytime.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

interface_primme.o: interface_primme.c interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

critical.o: critical.c critical.h geometry.h prob.h residual.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

plot.o: plot.c plot.h geometry.h prob.h residual.h critical.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

euler.o: euler.c euler.h geometry.h prob.h residual.h critical.h mytime.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

arpack.o: arpack.c arpack.h geometry.h prob.h residual.h mytime.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

%.o: %.c %.h
	$(CC) $(COPT) -c $< -o $@

# Verifie la matrice generee contre les fichiers CSR de reference du prof.
check: main
	./main -m 8 --check --no-solve

clean:
	rm -f *.o main main.exe

.PHONY: default check clean
