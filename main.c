#include "aux.h"
#include "householder.h"
#include "trelica.h"
#include <stdio.h>

void teste_matriz_ortogonal(matrix* V, int n){
    matrix* test = zeros(n);

    //transpoe V e guarda em test
    transpose(test, V, n);

    multiply_sq_matrix(test, test, V, n);
    int equals = 1;
    for (int i = 0; i < n ; i++)
        for (int j = 0; j < n; j++)
            if (i == j) {
                if (fabs(test->elem[i][j] - 1.0) > 1e-5){ //ignorar pequenas diferenças
                    equals=0;         
                }
            }else {
                if (fabs(test->elem[i][j] - 0.0) > 1e-5) //ignorar pequenas diferenças
                    equals=0;
            }
    if (equals == 1) 
        printf("    >A matriz dos autovetores é ortogonal.\n");
    else
        printf("    >A matriz dos autovetores não é ortogonal! Pânico!!!\n");
}

void teste_autovec(double* autovec, matrix* A, double autoval){
    int n = A->rows;
    double* y = calloc (n, sizeof(double));

    multiply_matrix_vector(y, A, autovec);

    int equals = 1;
    for (int i=0; i <n ; i++){
        y[i] /= autovec[i];
        if(fabs(y[i]-autoval) >= 1e-5)
            if ((fabs(autovec[i]) >= 1e-5) && (fabs(y[i]) >= 1e-5)) {
                equals = 0;            
            } //ignora divisão de zero por zero
    }

    if (equals == 1)
        printf("    >O autovetor calculado correspondente ao autovalor %3.2lf calculado é autovetor da matriz\n", autoval);
    else
        printf("    >O autovetor calculado do autovalor %3.2f não é autovetor da matriz A! Pânico!\n", autoval);
}

void teste(matrix* A){
    int n = A->rows;
    matrix* H = zeros(n);
    matrix* T = zeros(n);
    householder(A, T, H);

    printf(">Matriz A:");
    print_matrix(A, n);
    printf("\n>Matriz T");
    print_matrix(T, n);
    printf("\n>Matriz H:");
    print_matrix(H, n);

    matrix* V = zeros(n);
    matrix* AV = zeros(n);
    QRalgorithm(T, H, V, AV, n);

    printf("\n>Autovalores:\n");
    for (int i = 0; i<n; i++) 
        printf("%3.2f, ", AV->elem[i][i]);

    printf("\n>Autovetores:");
    print_matrix(V, n);    

    printf(">Testando se a matriz de autovetores é ortogonal\n");
    teste_matriz_ortogonal(V, n);
    //teste
    

    printf(">Testando se os autovetores multiplicando a matriz dão um multiplo do proprio autovetor\n");
    double*autovec = calloc(n, sizeof(double));
    for (int j=0; j < n; j++){
        for (int i = 0 ; i<n; i++){
            autovec[i] = V->elem[i][j];    
        }    
        teste_autovec(autovec, A, AV->elem[j][j]);
    }
    return;
}


int main(){
    printf(">>>EP1 MAP3121 - 2021<<<\n>Rodrigo Dias Siqueira (9425595) e Vitor Angelo Nunes(9836822)\n");
    int run = 1;
        while (run == 1) {
        printf("\n\n>Insira o teste que deseja realizar (a, b, c ou 0 para fechar)");
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

            teste(A);            
            fclose(f);
                   
        } else if (item == 'b') {
            FILE* f = fopen("input-b", "r");
            int n;
            fscanf(f, "%d", &n);

            matrix* A = zeros(n);

            for (int i = 0; i < n; i++)
                for (int j = 0; j<n; j++) 
                    fscanf(f, "%lf", &A->elem[i][j]);

            teste(A);            
            fclose(f);

        } else if (item =='c') {
            item_c();            
        } else if (item == '0')
            run = 0;
    }
    return 1;
}
