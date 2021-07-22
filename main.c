#include "householder.h"
#include <stdio.h>

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
    transpose(H, H, n);
    QRalgorithm(T, H, V, AV, n, 0);

    printf("\nAutovalores:\n");
    for (int i = 0; i<n; i++) 
        printf("%3.2f, ", AV->elem[i][i]);

    printf("\nAutovetores:");
    print_matrix(V, n);    
    return;
}



int main(){

    FILE* f = fopen("input-test", "r");
    int n;
    fscanf(f, "%d", &n);

    matrix* A = zeros(n);

    for (int i = 0; i < n; i++)
        for (int j = 0; j<n; j++) 
            fscanf(f, "%lf", &A->elem[i][j]);

    testecuca(A);
    fclose(f);
    return 1;

}