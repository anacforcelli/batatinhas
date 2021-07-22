#include "householder.h"
#include "aux.h"
#include <math.h>
#include <stdlib.h>

void householder (matrix* A){ 

    //Algoritmo baseado na seção 9.4 da 10a edição do Numerical Analysis, BURDEN-FAIRES.
    int n = A->rows;

    double*A_vec = malloc(((pow(n, 2)+n)/2)*sizeof(double));
    symm_matrix_to_vector(A_vec, A); //A_vec é o vetor q representa a matriz A,
                                     //coluna a coluna a partir da diagonal principal

    double*P = malloc(((pow(n, 2)+n)/2)*sizeof(double)); //P é o vetor que reresenta a matriz
                                                         //das transformacoes

    //step 1
    for (int k = 0; k < n-2; k++) { 
        
        //step 2 -check
        double q = 0;
        for (int j = k+1; j < n ; j++) {  
            q += pow(A_vec[convert_indices(j, k, n)], 2);
        }

        //step 3 -check
        double alpha;
        if (A_vec[convert_indices(k+1, k, n)] == 0) { 
            alpha = -sqrt(q);
        }else{
            alpha = -(sqrt(q) * A_vec[convert_indices(k+1, k, n)]) / fabs(A_vec[convert_indices(k+1, k, n)]);
        }

        //step 4 - check
        double RSQ = pow(alpha, 2) - alpha*A_vec[convert_indices(k+1, k, n)];
        
        //step 5 - check
        double* v = malloc((n-k)*sizeof(double));
        v[k]= 0;
        v[k+1] = A_vec[convert_indices(k+1, k, n)] - alpha;
        for (int j = k+2; j < n; j++) {
            v[j]=A_vec[convert_indices(j, k, n)];
        }

        double* w = malloc((n-k)*sizeof(double)+1);
        for (int i = 0; i < n-k ; i++)
            w[i]=(1/sqrt(2*RSQ))*v[i];
          
        //step 6 - check
        double* u = malloc((n-k)*sizeof(double)+1);
        for (int j = k; j < n; j++) {
            u[j]=0;
            for (int i = k+1 ; i < n ; i++) {
                u[j]+= A_vec[convert_indices(j, i, n)] * v[i];
            }
            u[j]/=RSQ;
        }

        //step 7
        double PROD = 0;
        for (int i = k+1; i < n; i++) {
            PROD += v[i]*u[i];
        }

        //step 8
        double* z = malloc((n-k)*sizeof(double)+1);
        for (int j = k; j < n; j++) {
            z[j] = u[j] - (0.5 * PROD/RSQ)*v[j];
        }

        //step 9   
        for (int l = k+1; l < n-1; l++){

            //step 11
            A_vec[convert_indices(l, l, n)] -= 2*v[l]*z[l]; 

            //step 10
            for (int j = l+1; j < n; j++) {
                A_vec[convert_indices(j, l, n)] -= v[l]*z[j] + v[j]*z[l];                
            } 
        }
        //step 12
        A_vec[convert_indices(n-1, n-1, n)] -= 2*v[n-1]*z[n-1];
        
        //step 13
        for (int j = k+2; j < n; j++){
            A_vec[convert_indices(k, j, n)] = 0;
        }

        //step 14
        A_vec[convert_indices(k+1, k, n)] -= v[k+1]*z[k];

        free(z);
        free(u);
        free(v);
        free(w);
    }
}