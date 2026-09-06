# librairies de PRIMME
LIBP = -L./primme/ -lprimme
# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ 
# toutes les librairies
LIB = $(LIBP) -lm -lblas -llapack

COPT = -O3 -Wall

default: main

clean: 
	rm *.o 
	rm main

main: main.o prob.o time.o interface_primme.o
	cc $(COPT) $^ -o $@ $(LIB)

main.o: main.c prob.h time.h interface_primme.h
	cc $(COPT) -c $< -o $@ $(INCP)

%.o: %.c %.h
	cc $(COPT) -c $< -o $@ $(INCP)


