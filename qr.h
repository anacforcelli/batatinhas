#include "aux.h"

#pragma once

void givens(matrix* A, matrix* Q, matrix* R, int n);
int QRalgorithm(matrix* A, matrix*HT ,matrix* V, matrix* AV, int n);