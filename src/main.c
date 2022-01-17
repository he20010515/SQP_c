#include "matrix.h"
#include "linearequations.h"
int main(void)
{
    Matrix *A = matrix_alloc(3, 3);
    Matrix *B = matrix_alloc(3, 3);

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            A->matrix_entry[i][j] = i + j;
            B->matrix_entry[i][j] = i * j;
        }
    }

    matrix_fill(A);
    B = matrix_inv(A);
    matrix_print(B);
}