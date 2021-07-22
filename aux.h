#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#pragma once

struct matrix {
  int rows, cols;
  double **elem ;
};


typedef struct matrix matrix;

void* meu_calloc(int times, int size);

matrix* zeros(int n);

void identity(matrix *res, int n);

void transpose(matrix*res, matrix* input, int n);

void copy_matrix(matrix*res, matrix*input, int n);
void multiply_matrix_scalar(matrix*out, matrix* in, double scalar, int n);

void multiply_sq_matrix(matrix*res, matrix*A, matrix*B, int n);

void print_matrix(matrix* A, int n);

void symm_matrix_to_vector(double* v, matrix* A);

int convert_indices(int i,int j, int n);

void vector_to_symm_matrix(double*v, matrix* A);

