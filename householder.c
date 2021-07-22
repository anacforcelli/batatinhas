#include "householder.h"
#include <math.h>
#include <stdlib.h>

void householder (matrix* A){ 

    //Algoritmo baseado na seção 9.4 da 10a edição do Numerical Analysis, BURDEN-FAIRES.
    //Adaptado para a representação de A por um vetor
    int n = A->rows;
    double*A_vec = malloc(((pow(n, 2)+n)/2)*sizeof(double)); ;
    symm_matrix_to_vector(A_vec, A); //A_vec é o vetor q representa a matriz A,
                                     //coluna a coluna a partir da diagonal principal

    int v_iter = 0; //indice do primeiro elemento da submatriz da iteracao k no vetor

    //step 1
    for (int k = 0; k < n-2; k++) { 
        
        //step 2 -check
        double q = 0;
        for (int j = 1; j < n-k ; j++) {  
            q += pow(A_vec[v_iter+j], 2);
        }

        //step 3 -check
        double alpha;
        if (A_vec[v_iter+1] == 0) { 
            alpha = -sqrt(q);
        }else{
            alpha = -(sqrt(q) * A_vec[v_iter+1]) / fabs(A_vec[v_iter+1]);
        }

        //step 4 - check
        double RSQ = pow(alpha, 2) - alpha*A_vec[v_iter+1];
        
        //step 5 - check
        double* v = malloc((n-k)*sizeof(double));
        v[k]= 0;
        v[k+1] = A_vec[v_iter+1] - alpha;
        for (int j = 2; j < n-k  ; j++) {
            v[j]=A_vec[v_iter+j];
        }

        double* w = malloc((n-k)*sizeof(double));
        for (int i = 0; i < n-k ; i++)
            w[i]=(1/sqrt(2*RSQ))*v[i];

        //step 6
        double* u = malloc((n-k)*sizeof(double));
        int v_iter_step_6 = v_iter;
        for (int j = 0; j<n-k; j++) {
            u[j]=0;
            for (int i = 1 ; i < n-k-j ; i++) {
                v_iter_step_6++;
                u[j]+= A_vec[v_iter_step_6] * v[i+k];
            }
            u[j]/=RSQ;
            v_iter_step_6++;
        }

        //step 7
        double PROD = 0;
        for (int i = k+1; i<n; i++) {
            PROD += v[i]*u[i];
        }

        //step 8
        double*z = malloc((n-k)*sizeof(double));
        for (int j=k; j<n; j++) {
            z[j] = u[j] - (0.5 * PROD/RSQ)*v[j];
        }

        //step 9
        int v_iter_step_9 = v_iter + (n-k) +1; //1 elemento após o elemento da diagonal principal     
        for (int l = k+1; l < n-1; l++){

            //step 11
            A_vec[v_iter_step_9]-=2*v[l]*z[l]; 

            //step 10
            for (int j = l+1; j<n; j++) {
                A_vec[v_iter_step_9] -= v[l]*z[j] + v[j]*z[l];
                v_iter_step_9++;
            }
            v_iter_step_9++;       
        }
        A_vec[v_iter_step_9] -= 2*v[n-k-1]*z[n-k-1];

        A_vec[v_iter+1] -= v[k+1]*z[k];
        //atualizacao do iterador 
        v_iter += n-k;
   }

}