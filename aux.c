#include "aux.h"

matrix* zeros(int n){
    matrix *m;
    int i;
    m = malloc(sizeof(matrix));
    m->rows = n;
    m->cols = n;
    m->elem = malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        m->elem[i] = calloc(n, sizeof(double));
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