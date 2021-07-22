#include "aux.h"

void givens(matrix* A, matrix* Q, matrix* R, int n){
    matrix* Q_aux = zeros(n);
    identity(Q_aux, n);
    matrix* R_aux = zeros(n);
    copy_matrix(R_aux, A, n);

    for (int i = 0; i < n-1; i++) {
        matrix* G = zeros(n);
        identity(G, n);
        double r, s, c;

        if (fabs(R_aux->elem[i+1][i])<=fabs(R_aux->elem[i][i])){
			r = -(R_aux->elem[i+1][i]/R_aux->elem[i][i]);
			c = 1.0/hypot(1.0, r);
			s = c*r;
        } else {
			r = -(R_aux->elem[i][i]/R_aux->elem[i+1][i]);
			s = 1.0/hypot(1.0, r);
			c = s*r;
        }
		G->elem[i][i] = c;
		G->elem[i+1][i+1] = c;
		G->elem[i][i+1] = -s;
		G->elem[i+1][i] = s;

		multiply_sq_matrix(R_aux, G, R_aux, n); //mesmo que A = GA

        transpose(G, G ,n);

		multiply_sq_matrix(Q_aux, Q_aux, G, n); // Q = Q * Gt

    }
    copy_matrix(Q, Q_aux, n);
    copy_matrix(R, R_aux, n);
}


int QRalgorithm(matrix* A, matrix* V, matrix* AV, int n, int shift){
    matrix* A_QR = zeros(n);
    copy_matrix(A_QR, A, n);
    matrix* V_aux = zeros(n);
	identity(V_aux, n);

	double eps = 0.000001;

    int	k = 0;

	for (int m = n-1; m > 0; m--) {
		double mu = 0;

		while (fabs(A_QR->elem[m][m-1]) > eps){
            matrix* Q = zeros(n);
            matrix* R = zeros(n);            

            if (k > 0 && shift==1){
                double d = 0.5 * (A_QR->elem[m-1][m-1] - A_QR->elem[m][m]);
                if (d>=0)
                    mu = A_QR->elem[m][m] + d - hypot(d, A_QR->elem[m][m-1]);
                else
                    mu = A_QR->elem[m][m] + d + hypot(d, A_QR->elem[m][m-1]);                     
            }
            
            for (int i=0; i<n; i++)
                A_QR->elem[i][i] -= mu;

            givens(A_QR, Q, R, n);
			multiply_sq_matrix(A_QR, R, Q, n);

            for (int i=0; i<n; i++)
                A_QR->elem[i][i] += mu;

			multiply_sq_matrix(V_aux, V_aux, Q, n);
			k++;
        }
        A_QR->elem[m][m-1] = 0.0;
    }
    copy_matrix(V, V_aux, n);
    copy_matrix(AV, A_QR, n);
    return k;
}
