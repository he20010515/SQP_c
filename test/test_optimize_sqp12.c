/*
 * @Author: HeYuwei
 * @Date: 2022-04-26 18:44:11
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-06-29 12:24:21
 * @FilePath: \SQP_c\test\test_optimize_sqp4.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "sqp.h"
#include "math.h"
#include "elog.h"
double _fun(Vector *X)
{
    double x1 = X->entry[0];
    double x2 = X->entry[1];
    double x3 = X->entry[2];
    double x4 = X->entry[3];

    return exp(x1 * x2) +
           5 * x4 -
           x1 * x3 * x4 -
           x1 * x2 * x2 * x3 * x3 * x3 +
           x2 * x2;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    double x3 = x->entry[2];
    double x4 = x->entry[3];

    y->entry[0] = 15 - x1;
    y->entry[1] = x2 - 4;
    y->entry[2] = 7 - x2;
    y->entry[3] = x1 * x2 + x3 - x4 - 6;
    y->entry[4] = 12 - x1 * x2 - x3 + x4;
    y->entry[5] = x3 - 5;
    y->entry[6] = 8 - x3;
    y->entry[7] = x4;
    y->entry[8] = 6 - x4;
    return;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 4, 9);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(4, 9, 0, 9, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 4);
    Vector *x0 = vector_alloc(4);
    Vector *xstar = vector_alloc(4);
    Vector *lambda0 = vector_alloc(9);

    vector_fill_const(x0, 1.0);

    vector_fill_const(lambda0, 0.0);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);
    printf(" min f = %lf\n", _fun(xstar));
    return 0;
}
