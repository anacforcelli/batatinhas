#include "householder.h"

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
   
    householder(A);
    return 1;
}