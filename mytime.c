#include <time.h>

/* retourne le temps CPU depuis le début de l'exécution du programme */
double mytimer_cpu(void){
    return (double) clock() / CLOCKS_PER_SEC;
}

#include <sys/time.h> 
/* retourne le temps d'horloge depuis le début de l’exécution du programme */
double mytimer_wall(){
        struct timeval dummy;
        gettimeofday( &dummy, NULL );
        return (double) ( (double) dummy.tv_sec +  ((double) dummy.tv_usec)/1000000.0 );
}

