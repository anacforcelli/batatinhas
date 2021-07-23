#include "aux.h"
#include <stdio.h>

void item_c (){
    FILE*  f = fopen("input-c", "r");
    int nos, nos_n_fixos, barras;
    fscanf(f, "%d %d %d", &nos, &nos_n_fixos, &barras);
    int ro, E;
    double A;
    fscanf(f, "%d %f %d", &ro, &A, &E);
}