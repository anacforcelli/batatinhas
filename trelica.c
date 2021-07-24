#include "aux.h"
#include "trelica.h"
#include "householder.h"
#include <stdio.h>
#include <stdlib.h>

void item_c (){

    //montar a estrutura da trelica a partir do arquivo de texto
    FILE*  f = fopen("input-c", "r");
    trelica* t = malloc(sizeof(trelica));
    
    fscanf(f, "%d %d %d", &t->n_nos, &t->n_nos_nao_fixos, &t->n_barras);

    t->barras = calloc_wrapper(t->n_barras,sizeof(barra));
    for (int i = 0; i < t->n_barras; i++)
        t->barras[i] = calloc_wrapper(1, sizeof(barra));
    t->nos = calloc_wrapper(t->n_nos,sizeof(no));
    for (int i = 0; i < t->n_nos; i++)
        t->nos[i] = calloc_wrapper(1, sizeof(no));
    
    fscanf(f, "%d %lf %d", &t->ro, &t->A, &t->E);
    for (int i =0 ; i < t->n_barras; i++) {
        int no_1, no_2;
        fscanf(f, "%d %d %lf %lf", 
        &no_1,
        &no_2,
        &t->barras[i]->angulo,
        &t->barras[i]->comprimento);
        t->barras[i]->no_1 = t->nos[no_1-1];
        t->barras[i]->id_no_1 = no_1;
        t->barras[i]->no_2 = t->nos[no_2-1];
        t->barras[i]->id_no_2 = no_2;
    }
    
    //inicializar as massas atribuidas a cada no
    for (int i=0; i<t->n_nos; i++)
        t->nos[i]->massa=0;

    for (int i=0 ; i < t->n_barras; i++)
        t->barras[i]->no_1->massa  +=  0.5 * t->ro * t->A * t->barras[i]->comprimento;
    

    //matriz diagonal M = M^(-1/2) 
    double*M = calloc(2*t->n_nos_nao_fixos,sizeof(double));
    for( int i = 0; i < t->n_nos_nao_fixos; i++)
        M[2*i] = pow(t->nos[i]->massa, -0.5);
    

    //inicializar matriz Kij de cada barra (está muito feio, mas nao consegui pensar de outro jeito)
    for (int k = 0; k < t->n_barras; k++) {
        t->barras[k]->Kij = zeros(4);
        double L = t->barras[k]->comprimento;
        double theta = t->barras[k]->angulo;
        t->barras[k]->Kij->elem[0][0] = ((t->A*t->E)/L) * cos(theta)*cos(theta);
        t->barras[k]->Kij->elem[0][1] = ((t->A*t->E)/L) * cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[0][2] = ((t->A*t->E)/L) * -cos(theta)*cos(theta);
        t->barras[k]->Kij->elem[0][3] = ((t->A*t->E)/L) * -cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[1][0] = ((t->A*t->E)/L) * cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[1][1] = ((t->A*t->E)/L) * sin(theta)*sin(theta);
        t->barras[k]->Kij->elem[1][2] = ((t->A*t->E)/L) * -cos(theta)*cos(theta);
        t->barras[k]->Kij->elem[1][3] = ((t->A*t->E)/L) * -sin(theta)*sin(theta);
        t->barras[k]->Kij->elem[2][0] = ((t->A*t->E)/L) * -cos(theta)*cos(theta);
        t->barras[k]->Kij->elem[2][1] = ((t->A*t->E)/L) * -cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[2][2] = ((t->A*t->E)/L) * cos(theta)*cos(theta);
        t->barras[k]->Kij->elem[2][3] = ((t->A*t->E)/L) * cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[3][0] = ((t->A*t->E)/L) * -cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[3][1] = ((t->A*t->E)/L) * -sin(theta)*sin(theta);
        t->barras[k]->Kij->elem[3][2] = ((t->A*t->E)/L) * cos(theta)*sin(theta);
        t->barras[k]->Kij->elem[3][3] = ((t->A*t->E)/L) * cos(theta)*cos(theta);
    }
    
    //inicializar matriz K de rigidez total da treliça
    t->K = zeros(2*t->n_nos_nao_fixos);
    for (int k = 0; k<t->n_nos_nao_fixos; k++) {
        int i = t->barras[k]->id_no_1;
        int j = t->barras[k]->id_no_2;
        if ((i != 10 && j != 13)||(i!=11&&j!=12)) {//se nao forem as barras 11,14 ou 12, 13
            t->K->elem[2*i-1][2*i-1] = t->barras[k]->Kij->elem[0][0];
            t->K->elem[2*i-1][2*i]   = t->barras[k]->Kij->elem[0][1];
            t->K->elem[2*i-1][2*j-1] = t->barras[k]->Kij->elem[0][2];
            t->K->elem[2*i-1][2*j]   = t->barras[k]->Kij->elem[0][3];
            t->K->elem[2*i][2*i-1]   = t->barras[k]->Kij->elem[1][0];
            t->K->elem[2*i][2*i]     = t->barras[k]->Kij->elem[1][1];
            t->K->elem[2*i][2*j-1]   = t->barras[k]->Kij->elem[1][2];
            t->K->elem[2*i][2*j]     = t->barras[k]->Kij->elem[1][3];
            t->K->elem[2*j-1][2*i-1] = t->barras[k]->Kij->elem[2][0];
            t->K->elem[2*j-1][2*i]   = t->barras[k]->Kij->elem[2][1];
            t->K->elem[2*j-1][2*j-1] = t->barras[k]->Kij->elem[2][2];
            t->K->elem[2*j-1][2*j]   = t->barras[k]->Kij->elem[2][3];
            t->K->elem[2*j][2*i-1]   = t->barras[k]->Kij->elem[3][0];
            t->K->elem[2*j][2*i]     = t->barras[k]->Kij->elem[3][1];
            t->K->elem[2*j][2*j-1]   = t->barras[k]->Kij->elem[3][2];
            t->K->elem[2*j][2*j]     = t->barras[k]->Kij->elem[3][3];
        } /*else {
            t->K->elem[2*i-1][2*i-1] = t->barras[k]->Kij->elem[0][0];
            t->K->elem[2*i-1][2*i]   = t->barras[k]->Kij->elem[0][1];
            t->K->elem[2*i-1][2*j-1] = t->barras[k]->Kij->elem[0][2];
            t->K->elem[2*i-1][2*j]   = t->barras[k]->Kij->elem[0][3];
            t->K->elem[2*i][2*i-1]   = t->barras[k]->Kij->elem[1][0];
            t->K->elem[2*i][2*i]     = t->barras[k]->Kij->elem[1][1];
            t->K->elem[2*i][2*j-1]   = t->barras[k]->Kij->elem[1][2];
            t->K->elem[2*i][2*j]     = t->barras[k]->Kij->elem[1][3];
            t->K->elem[2*j-1][2*i-1] = t->barras[k]->Kij->elem[2][0];
            t->K->elem[2*j-1][2*i]   = t->barras[k]->Kij->elem[2][1];
            t->K->elem[2*j-1][2*j-1] = t->barras[k]->Kij->elem[2][2];
            t->K->elem[2*j-1][2*j]   = t->barras[k]->Kij->elem[2][3];
            t->K->elem[2*j][2*i-1]   = t->barras[k]->Kij->elem[3][0];
            t->K->elem[2*j][2*i]     = t->barras[k]->Kij->elem[3][1];
            t->K->elem[2*j][2*j-1]   = t->barras[k]->Kij->elem[3][2];
            t->K->elem[2*j][2*j]     = t->barras[k]->Kij->elem[3][3];
        }*/

        //multiplicar k por M^-0.5 pela esquerda e direita
        for(int i =0; i < t->n_nos_nao_fixos; i++)
            for (int j = 0; j < t->n_nos_nao_fixos; j++) 
                t->K->elem[i][j] *= M[i]*M[j];

        matrix* T = zeros(2*t->n_nos_nao_fixos);
        matrix* H = zeros(2*t->n_nos_nao_fixos);
        
        householder(t->K, T, H);

        matrix*V = zeros(2*t->n_nos_nao_fixos);
        matrix*AV = zeros(2*t->n_nos_nao_fixos);
        QRalgorithm(T, H,V,AV, 2*t->n_nos_nao_fixos);

        printf("Os modos de vibração de maior frequencia sao:");

        for (int i = 0; i < 5; i++){
            printf("%3.2lf", sqrt(AV->elem[i][i]));
        }
        
    }
}
    


