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

# toutes les librairies
LIB = $(LIBP) $(LIBBLAS) -lm

COPT = -O3 -Wall -Wextra

OBJ = main.o geometry.o prob.o residual.o mytime.o interface_primme.o

default: main

main: $(OBJ)
	$(CC) $(COPT) $^ -o $@ $(LIB)

main.o: main.c main.h geometry.h prob.h residual.h mytime.h interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

interface_primme.o: interface_primme.c interface_primme.h
	$(CC) $(COPT) -c $< -o $@ $(INCP)

%.o: %.c %.h
	$(CC) $(COPT) -c $< -o $@

# Verifie la matrice generee contre les fichiers CSR de reference du prof.
check: main
	./main -m 8 --check --no-solve

clean:
	rm -f *.o main main.exe

.PHONY: default check clean
