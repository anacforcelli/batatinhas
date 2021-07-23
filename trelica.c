#include "aux.h"
#include "trelica.h"
#include <stdio.h>
#include <stdlib.h>

void item_c (){

    //montar a estrutura da trelica a partir do arquivo de texto
    FILE*  f = fopen("input-c", "r");
    trelica* t = malloc(sizeof(trelica));
    
    fscanf(f, "%d %d %d", &t->n_nos, &t->n_nos_nao_fixos, &t->n_barras);

    t->barras = malloc(t->n_barras*sizeof(barra));
    t->nos = malloc(t->n_nos*sizeof(no));
    
    fscanf(f, "%d %lf %d", &t->ro, &t->A, &t->E);
    for (int i =0 ; i < t->n_barras; i++) {
        int no_1, no_2;
        fscanf(f, "%d %d %lf %lf", 
        &no_1,
        &no_2,
        &t->barras[i]->angulo,
        &t->barras[i]->comprimento);
        t->barras[i]->no_1 = t->nos[no_1-1];
        t->barras[i]->no_2 = t->nos[no_2-1];
    }
    
    //inicializar as massas atribuidas a cada no
    for (int n=0; n<t->n_nos; n++){
        t->nos[n]->massa=0;
    }
    
    for (int i=0 ; i < t->n_barras; i++){
        t->barras[i]->no_1->massa  +=  0.5 * t->ro * t->A * t->barras[i]->comprimento;
    }

    



}