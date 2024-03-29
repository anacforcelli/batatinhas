#include "aux.h"
#include <stdio.h>
#include <stdlib.h>


void* calloc_wrapper(int times, int size) { 
    //printf("calloc(%d)\n", times*size); 
    void* a = malloc(times*size);
    if(!a){
        printf("MALLOC FAILED");
        exit(1);
    }
    return a;
}

matrix* zeros(int n){
    matrix *m;
    int i;
    m = calloc_wrapper(1, sizeof(matrix));
    m->rows = n;
    m->cols = n;
    m->elem = calloc_wrapper(n, sizeof(double *));
    for (i = 0; i < n; i++) {
        m->elem[i] = calloc_wrapper(n, sizeof(double));
    }
    return m;    
}

void identity(matrix*res, int n){
    for (int i = 0; i < n; i++){   
        for(int j = 0; j < n; j++){
            if (j != i)
                res->elem[i][j] = 0;
            else
                res->elem[i][j] = 1;
        }
    }   
}

void transpose(matrix*res, matrix* input, int n){
    matrix* aux = zeros(n);
    for (int i = 0; i < n; i++){   
        for(int j = 0; j < n; j++){
            aux->elem[i][j] = input->elem[j][i];
        }
    }   
    copy_matrix(res, aux, n);
}

void copy_matrix(matrix*res, matrix*input, int n){    
    for (int i = 0; i < n; i++){   
        for(int j = 0; j < n; j++){
            res->elem[i][j] = input->elem[i][j];
        }
    }   
}

void multiply_matrix_scalar(matrix*out, matrix* in, double scalar, int n){
    for (int i = 0; i < n; i++){   
        for(int j = 0; j < n; j++){
            out->elem[i][j] = in->elem[i][j]*scalar;
        }
    }   
}

void multiply_sq_matrix(matrix*res, matrix* A, matrix*B, int n){   
    matrix*aux = zeros(n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            for (int k = 0; k < n; k++){
                aux->elem[i][j] += A->elem[i][k]*B->elem[k][j];
            }
        }
    }   
    copy_matrix(res, aux, n);
}

void multiply_matrix_vector(double* res, matrix* A, double* v){
    int n = A->rows;

    for (int i = 0; i<n; i++){
        res[i] = 0;
        for (int j=0; j<n; j++)
            res[i] += A->elem[i][j]*v[j];
    }
}


void print_matrix(matrix* A, int n){
    for (int i = 0; i < n; i++){   
        printf("\n[");
        for(int j = 0; j < n; j++){
            printf("%2.2f, ", A->elem[i][j]);
        }
        printf("]\n");
    }   
}

//lê a matriz triangular superior, armazenando-a linha
//a linha a partir do elemento da diagonal principal.
void symm_matrix_to_vector(double* A_vec, matrix* A){
    int k = 0;
    int n = A->rows;
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n ; j++) {
            A_vec[k]=A->elem[i][j];
            k++;
        }
    }
}

//converte os indices de uma matriz comum para acessar o vetor
int convert_indices(int i,int j, int n){
    if (i>=n || i<0 || j>=n || j<0 ){
        printf("errouuuuu");
        return -1;
    }
    if (i > j){
        int temp = i;
        i = j;
        j = temp;
    }
    int x = i*n + j - (1+i)*(i)/2;
    return x;
}
   
void vector_to_symm_matrix(double*v, matrix* A, int n){
    for(int i=0; i<n; i++)
        for (int j=0 ; j<n ; j++) {
           A->elem[i][j] = v[convert_indices(i, j, n)];
           A->elem[j][i] = v[convert_indices(j, i, n)];
        }     
}
