#include "matrix.h"
#include "vector.h"

int optimize_qp(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy);