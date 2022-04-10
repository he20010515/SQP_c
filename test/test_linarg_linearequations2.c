#include "matrix.h"
#include "vector.h"
#include "linear_equations.h"

int main(int argc, char const *argv[])
{
    const int N1 = 3;
    double A_2darray[3][3] = {
        {8., -3., 2.},
        {4., 11., -1.},
        {6., 3., 12.},
    };
    Matrix *A = matrix_alloc(N1, N1);
    array_2_matrix((double *)A_2darray, N1, N1, A);
    Vector *b = vector_alloc(N1);
    b->entry[0] = 20;
    b->entry[1] = 33;
    b->entry[2] = 36;

    Vector *x = vector_alloc(N1);
    Vector *x0 = vector_alloc(N1);
    x0->entry[0] = 0;
    x0->entry[1] = 0;
    x0->entry[2] = 0;

    linear_equation_jacobi(A, b, 0.0001, x0, x);
    Vector *x2 = vector_alloc(N1);
    vector_print(x);
    linear_equation_gauss_sidel(A, b, 0.0001, x0, x2);
    vector_print(x2);
    Vector *x3 = vector_alloc(N1);
    linear_equation_sor(A, b, 0.0001, x0, x3, 0.5);
    vector_print(x3);


    return 0;
}
