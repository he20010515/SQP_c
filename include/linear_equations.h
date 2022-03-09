#include "vector.h"
#include "matrix.h"

void linear_equation_jacobi(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x);