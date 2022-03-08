#include "function.h"
#include "matrix.h"
#include "math.h"
#include "function.h"
#include "vector.h"
int main(void)
{
    const int N1 = ;
    Matrix *A = matrix_alloc(N1, N1);
    for (size_t i = 0; i < N1; i++)
    {
        for (size_t j = 0; j < N1; j++)
        {
            A->matrix_entry[i][j] = (double)i * j + 3.;
        }
    }
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
