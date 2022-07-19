/*
 * @Author: HeYuwei
 * @Date: 2022-07-18 20:34:36
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-18 20:40:23
 * @FilePath: \SQP_c\test\test_simplex.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "sqp.h"
#include "simplex.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    double A[2][7] = {{1, 1, 1, -1, -1, -1, 0}, {21, 19, -43, -21, -19, 43, -1}};
    Matrix *A_mat = matrix_alloc(2, 7);
    array_2_matrix((double *)A, 2, 7, A_mat);

    Vector *c = vector_alloc(7);
    vector_fill_const(c, 1.0);

    Vector *b = vector_alloc(2);
    b->entry[0] = 2.30926e-013;
    b->entry[1] = 672.333;
    Vector *x = vector_alloc(7);

    _linprog_simplex(c, A_mat, b, 10, 1e-9, 0, x);
    vector_print(x);

    return 0;
}
