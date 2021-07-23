#include "householder.h"
#include "trelica.h"
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
    QRalgorithm(T, H, V, AV, n);

    printf("\nAutovalores:\n");
    for (int i = 0; i<n; i++) 
        printf("%3.2f, ", AV->elem[i][i]);

    printf("\nAutovetores:");
    print_matrix(V, n);    
    return;
}


int main(){
    printf("EP1 MAP3121\nRodrigo Dias Siqueira e Vitor Angelo Nunes\n");
    int run = 1;
        while (run == 1) {
        printf("\n\n--insira o teste que deseja realizar --a b c: (ESCREVA 0 PARA FECHAR)");
        char item;
        scanf("%c", &item);
        if (item == 'a'){
            FILE* f = fopen("input-a", "r");
            int n;
            fscanf(f, "%d", &n);

            matrix* A = zeros(n);

            for (int i = 0; i < n; i++)
                for (int j = 0; j<n; j++) 
                    fscanf(f, "%lf", &A->elem[i][j]);

            testecuca(A);            
            fclose(f);
                   
        } else if (item == 'b') {
            FILE* f = fopen("input-b", "r");
            int n;
            fscanf(f, "%d", &n);

            matrix* A = zeros(n);

            for (int i = 0; i < n; i++)
                for (int j = 0; j<n; j++) 
                    fscanf(f, "%lf", &A->elem[i][j]);

            testecuca(A);            
            fclose(f);

        } else if (item =='c') {
            item_c();            
        } else if (item == '0')
            run = 0;
    }
    return 1;
}
