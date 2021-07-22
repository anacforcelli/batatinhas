#include "householder.h"

void testecuca(matrix* A){
    int n = A->rows;
    matrix* H = zeros(n);
    matrix*T = zeros(n);
    householder(A, T, H);

    printf("Matriz A:");
    print_matrix(A, n);
    printf("\nMatriz T");
    print_matrix(T, n);
    printf("\nMatriz H:");
    print_matrix(H, n);

    matrix* V = zeros(n);
    matrix* AV = zeros(n);
    QRalgorithm_modded(T, H, V, AV, n, 1);

    printf("Autovalores:");
    print_matrix(AV, n);
}

int main(){


    matrix* A = zeros(4);

    A->elem[0][0] = 4;
    A->elem[1][0] = 1;
    A->elem[0][1] = 1;
    A->elem[2][0] = -2;
    A->elem[0][2] = -2;
    A->elem[3][0] = 2;
    A->elem[0][3] = 2;
    A->elem[1][1] = 2;
    A->elem[2][1] = 0;
    A->elem[1][2] = 0;
    A->elem[3][1] = 1;
    A->elem[1][3] = 1;
    A->elem[2][2] = 3;
    A->elem[2][3] = -2;
    A->elem[3][2] = -2;
    A->elem[3][3] = -1;


    
    return 1;
}