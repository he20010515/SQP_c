#include "matrix.h"
#include "linearequations.h"
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

    matrix_copy(A, B); // B = A
    printf("MatrixA:\n");
    matrix_print(A);
    printf("MatrixA^-1:\n");
    matrix_invert(B); // B = B^-1
    matrix_print(B);
    matrix_print(matrix_multiply(A, B));
    matrix_print(matrix_multiply(B, A));
}