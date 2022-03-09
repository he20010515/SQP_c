#include "matrix.h"
#include "vector.h"
#include "linear_equations.h"
int main(int argc, char const *argv[])
{
    const int N1 = 4;
    double A_2darray[4][4] = {
        {1.0, 7.0, 8.0, 1},
        {1.0, 4.0, 1.0, 1},
        {2.0, 1.0, 7.0, 1},
        {2.0, 2.0, 7.0, 1},
    };

    Matrix *Ainv = matrix_alloc(N1, N1);
    Vector *Ainvb = vector_alloc(N1);
    Vector *b = vector_alloc(4);
    b->entry[0] = 0;
    b->entry[1] = 1;
    b->entry[2] = 2;
    b->entry[3] = 3;
    Vector *x = vector_alloc(4);
    Vector *x0 = vector_alloc(4);
    x0->entry[0] = 7.42;
    x0->entry[1] = 1;
    x0->entry[2] = -0.5;
    x0->entry[3] = -9.85;

    Matrix *A = matrix_alloc(N1, N1);
    array_2_matrix((double *)A_2darray, N1, N1, A);
    printf("linerfunctions Ax=b\n A = \n");
    matrix_print(A);
    printf("b = \n");
    vector_print(b);
    printf("A^{-1} = \n");
    matrix_inverse(A, Ainv);
    matrix_print(Ainv);
    printf("x find by A^{-1}b:\n");
    matrix_mutiply_vector(Ainv, b, Ainvb);
    vector_print(Ainvb);
    Vector *AAinvb = vector_alloc(N1);
    matrix_mutiply_vector(A, Ainvb, AAinvb);
    printf("A *(A^{-1}b) = ");
    vector_print(AAinvb);
    printf("====jacobi method===\n");
    printf("x0:\n");

    vector_print(x0);
    linear_equation_jacobi(A, b, 0.01, x0, x);
    printf("x:");
    vector_print(x);
    return 0;
}
