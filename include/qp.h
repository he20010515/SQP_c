#pragma once\

#include "matrix.h"
#include "vector.h"
#include "lp.h"

int optimize_qp_linear_constraints(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy);
int optimize_qp_active_set(const Matrix *G, const Vector *c, const LinearConstraints *cons, const Vector *x0, Vector *x_star, Vector *lam);

