#include "matrix.h"
#include <memory.h>
char *test(void)
{
    printf("test function from linearquation");
}

Matrix *linear_solve(Matrix *A, Matrix *b)
{
    if (A->col_size != A->row_size) // 是否是矩阵
    {
        return -1;
    }
    if (A->row_size != b->row_size) //A b 是否匹配
    {
        return -1;
    }
    if (b->col_size != 1) // b 是否为列向量
    {
        return -1;
    }
    Matrix *x = matrix_alloc(b->row_size, 1);
    Matrix *temp = matrix_alloc(b->row_size, 1);

    matrix_free(temp);
    return x;
}