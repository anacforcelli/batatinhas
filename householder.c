#include "householder.h"
#include "aux.h"
#include <math.h>
#include <stdlib.h>

//A partir de uma matriz A, retorna sua forma tridiagonal e a matriz H das transformações. 
void householder (matrix*A, matrix*T, matrix*H){ 

    //Algoritmo baseado na seção 9.4 da 10a edição do Numerical Analysis, BURDEN-FAIRES.
    int n = A->rows;

    double*A_vec = calloc(((pow(n, 2)+n)/2),sizeof(double));
    symm_matrix_to_vector(A_vec, A); //A_vec é o vetor q representa a matriz A,
                                     //coluna a coluna a partir da diagonal principal

    //como H não é simétrica, não é tão custoso utilizá-la como uma matriz cheia
    identity(H, n);
    
    //Hwi, no entanto, é simétrica.
    double* H_wi = calloc(((pow(n, 2)+n)/2),sizeof(double));

    //step 1
    for (int k = 0; k < n-2; k++) { 
        
        for (int i=0; i<n; i++){ //inicializado H_wi com identidade
            for (int j=0; j<n; j++) 
                H_wi[convert_indices(i, j, n)]=0.0;  
            H_wi[convert_indices(i, i, n)]=1.0; 
        }

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
        //notar que para todos os malloc's, os elementos antes de k 
        //nao serao utilizados, logo nao preciso inicializá-los
        double* v = malloc(n*sizeof(double));
        v[k]= 0;
        v[k+1] = A_vec[convert_indices(k+1, k, n)] - alpha;
        for (int j = k+2; j < n; j++) {
            v[j]=A_vec[convert_indices(j, k, n)];
        }
          
        //step 6 - check
        double* u = malloc(n*sizeof(double));
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
        double* z = malloc(n*sizeof(double));
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

        // calcular w
        double* w = malloc(n*sizeof(double));
        for (int i = k; i < n ; i++)
            w[i]=(1/sqrt(2*RSQ))*v[i];

        // com w, atualizar a matriz H_k
        for (int i=k+1; i<n; i++)
            for (int j=k+1; j<=i; j++)
                H_wi[convert_indices(i, j, n)] -=  2*w[i]*w[j]; 

        //atualizar a matriz H a cada iteracao
        matrix* H_wi_m = zeros(n);
        vector_to_symm_matrix(H_wi, H_wi_m);
        multiply_sq_matrix(H, H, H_wi_m, n);
        
        free(v);
        free(w);
        free(H_wi);
    }    
    vector_to_symm_matrix(A_vec, T);
}