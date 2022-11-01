/*
 * @Author: heyuwei he20010515@163.com
 * @Date: 2022-11-01 20:53:00
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-01 23:47:15
 * @FilePath: /SQP_c/src/linarg/pinv.c
 * @Description:矩阵伪逆计算
 *
 */
#include "matrix.h"

void matrix_pinv(const Matrix *A, Matrix *A_plus, double tol)
{
    if (!(A_plus->row_size == A->col_size AND A_plus->col_size == A->row_size))
    {
        terminate("matrix_pinv: size dont fit");
    }

    Matrix *U = matrix_alloc(A->row_size, A->row_size),
           *W = matrix_alloc(A->row_size, A->col_size),
           *V_T = matrix_alloc(A->col_size, A->col_size);

    Matrix *U_T = matrix_alloc(A->row_size, A->row_size),
           *W_inv = matrix_alloc(A->col_size, A->row_size),
           *V = matrix_alloc(A->col_size, A->col_size);
    matrix_singluar(A, U, W, V_T);
    matrix_transpose(U, U_T);
    matrix_transpose(V_T, V);
    int minmn = 0;
    if (A->row_size < A->col_size)
        minmn = A->row_size;
    else
        minmn = A->col_size;
    matrix_fill_const(W_inv, 0);
    for (int i = 0; i < minmn; i++)
    {
        if (W->matrix_entry[i][i] > tol)
        {
            W_inv->matrix_entry[i][i] = 1 / W->matrix_entry[i][i];
        }
        else
        {
            W_inv->matrix_entry[i][i] = 0.0;
        }
    }

    Matrix *V_W_inv = matrix_multiply(V, W_inv);
    Matrix *V_W_inv_U_T = matrix_multiply(V_W_inv, U_T);
    matrix_copy(V_W_inv_U_T, A_plus);

    matrix_free(U);
    matrix_free(W);
    matrix_free(V_T);
    matrix_free(U_T);
    matrix_free(W_inv);
    matrix_free(V);
    matrix_free(V_W_inv);
    matrix_free(V_W_inv_U_T);
    return;
}