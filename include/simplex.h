#pragma once
#include "vector.h"
#include "matrix.h"

int _linprog_simplex(const Vector *c, const Matrix *A, const Vector *b, int maxiter, double tol, int bland, Vector *x);
