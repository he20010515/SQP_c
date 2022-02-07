#include "matrix.h"
#include "linearequations.h"
#include "f2c.h"
#include "slsqp.h"

#define M_alloc_variable(x, type, value)    \
    type *x = (type *)malloc(sizeof(type)); \
    *x = value

int main(void)
{
    Matrix *A = matrix_alloc(3, 3);
    Matrix *B = matrix_alloc(3, 3);
    Matrix *C = matrix_callalloc(3);
    A->matrix_entry[0][0] = 4.;
    A->matrix_entry[0][1] = 2.;
    A->matrix_entry[0][2] = 1.;
    A->matrix_entry[1][0] = 12.;
    A->matrix_entry[1][1] = 3.;
    A->matrix_entry[1][2] = 2.;
    A->matrix_entry[2][0] = 4.;
    A->matrix_entry[2][1] = 5.;
    A->matrix_entry[2][2] = 2.;

    integer a = 1;

    // matrix_copy(A, B); // B = A
    // printf("MatrixA:\n");
    // matrix_print(A);
    // printf("MatrixA^-1:\n");
    // matrix_invert(B); // B = B^-1
    // matrix_print(B);
    // matrix_print(matrix_multiply(A, B));
    // matrix_print(matrix_multiply(B, A));
    int state = 99;
    M_alloc_variable(m, integer, 0);
    M_alloc_variable(w, integer, 1);

    printf("state%d", *m);
    printf("state%d", *w);

}