/*
 * @Author: HeYuwei
 * @Date: 2022-07-15 12:26:15
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-15 12:34:35
 * @FilePath: \SQP_c\test\test_optimize_qp_active_set_fromsqp.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "qp.h"
#include "matrix.h"
#include "vector.h"
#include "index_set.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    Matrix *G = matrix_alloc(3, 3);
    double G_array[3][3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    array_2_matrix((double *)G_array, 3, 3, G);

    Vector *c = vector_alloc(3);
    Vector *b = vector_alloc(2);
    c->entry[0] = 0.;
    c->entry[1] = 0.;
    c->entry[2] = -32.;

    b->entry[0] = -2.;
    b->entry[1] = 15.;

    Matrix *A = matrix_alloc(2, 3);
    double A_array[2][3] = {
        {1, 1, 1},
        {-2, -4, -2},
    };
    array_2_matrix((double *)A_array, 2, 3, A);
    LinearConstraints *con = linearconstraints_alloc(3, 2, 1, 1, A, b);
    Vector *x0 = vector_alloc(3);
    Vector *lambda = vector_alloc(2);
    Vector *xstar = vector_alloc(3);
    x0->entry[0] = 1.2;
    x0->entry[1] = 0.8;
    optimize_qp_active_set(G, c, con, NULL, xstar, lambda);
    vector_print(xstar);
    return 0;
}
