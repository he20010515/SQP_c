#include "function.h"
#include "matrix.h"
#include "math.h"
#include "function.h"
#include "vector.h"

#define ADD(x, y)

int main(void)
{
    const int N1 = 4;
    Matrix *A = matrix_alloc(N1, N1);

    double A_2darray[4][4] = {
        {1.0, 7.0, 8.0, 1},
        {1.0, 4.0, 1.0, 1},
        {2.0, 1.0, 7.0, 1},
        {2.0, 2.0, 7.0, 1},
    };
    array_2_matrix((double *)A_2darray, N1, N1, A);
    Matrix *A_inv = matrix_alloc(N1, N1);
    matrix_inverse(A, A_inv);
    printf("Matrix A:\n");
    matrix_print(A);
    printf("Matrix A^-1:\n");
    matrix_print(A_inv);
    Matrix *A_A_1 = matrix_multiply(A, A_inv);

    printf("Matrix A*A^-1:\n");
    matrix_print(A_A_1);
}
