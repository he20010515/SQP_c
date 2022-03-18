#include "qp.h"
#include "matrix.h"
#include "vector.h"

int main(int argc, char const *argv[])
{
    // min \frac12 x^{T}Hx + c^{T}x
    // s.t. Ax = b
    const int N1 = 2;
    double A_2darray[2][2] = {
        {-1., 2.0},
        {0., 1.}};
    Matrix *A = matrix_alloc(N1, N1);
    array_2_matrix((double *)A_2darray, N1, N1, A);
    Matrix *H = matrix_callalloc(2);
    Vector *b = vector_alloc(N1);
    b->entry[0] = 0;
    b->entry[1] = 0;

    Vector *c = vector_alloc(N1);
    c->entry[0] = 1;
    c->entry[1] = 2;
    Vector *x_star = vector_alloc(N1);
    double y;
    optimize_qp_linear_constraints(H, c, A, b, x_star, &y);
    vector_print(x_star);
    return 0;
}