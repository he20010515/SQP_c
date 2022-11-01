/*
 * @Author: heyuwei he20010515@163.com
 * @Date: 2022-11-01 21:03:56
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-01 23:44:20
 * @FilePath: /SQP_c/test/test_linarg_pinv.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "matrix.h"
int main(int argc, char const *argv[])
{
    int m = 5;
    int n = 1;
    Matrix *A = matrix_alloc(m, n);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A->matrix_entry[i][j] = i * (j + 1) + 1;
        }
    }
    Matrix *U = matrix_alloc(m, m);
    Matrix *W = matrix_alloc(m, n);
    Matrix *VT = matrix_alloc(n, n);
    printf("Matrix A");
    matrix_print(A);
    matrix_singluar(A, U, W, VT);

    printf("Matrix U");
    matrix_print(A);
    printf("Matrix W");
    matrix_print(W);
    printf("Matrix VT");
    matrix_print(VT);

    Matrix *temp1 = matrix_multiply(U, W);
    Matrix *temp2 = matrix_multiply(temp1, VT);
    printf("Matrix U*W*VT");
    matrix_print(temp2);

    Matrix *A_plus = matrix_alloc(A->col_size, A->row_size);
    printf("Matrix A^+\n");
    matrix_pinv(A, A_plus, 1e-8);
    matrix_print(A_plus);
    return 0;
}
