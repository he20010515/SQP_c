#include "vector.h"
#include "matrix.h"
#pragma once

void linear_equation_jacobi(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x);
void linear_equation_gauss_sidel(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x);
void linear_equation_sor(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x, const double omega);
