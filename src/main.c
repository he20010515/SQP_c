#include "matrix.h"
#include "linearequations.h"
int main(void)
{
    Matrix *A = matrix_alloc(2, 2);
    Matrix *B = matrix_alloc(2, 2);

    A->matrix_entry[0][0] = 1;
    A->matrix_entry[0][1] = 0;
    A->matrix_entry[1][0] = 10;
    A->matrix_entry[1][1] = 2;

    matrix_copy(A, B); // B = A
    matrix_print(A);
    matrix_print(B);

    matrix_invert(B); // B = B^-1
    matrix_print(B);
    Matrix *C = matrix_multiply(A, B); // C = A*B = A*B^-1 = A*A^-1;
    matrix_print(C);
}