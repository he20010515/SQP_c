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
    double x5 = X->entry[4];

    return pow(x1, 6) +
           7 * x1 * x2 * x2 +
           5 * pow(x3, 3.0) +
           6 * pow(x4, 4) -
           3 * pow(x5, 6) -
           4 * x3 * x4 -
           6 * x3 * x3 * x5 -
           x1 * x4;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    double x3 = x->entry[2];
    double x4 = x->entry[3];
    double x5 = x->entry[4];

    y->entry[0] = x1 - 3;
    y->entry[1] = 6 - x1;
    y->entry[2] = x2;
    y->entry[3] = 3 - x2;
    y->entry[4] = x1 * x2 - 4;
    y->entry[5] = 7 - x1 * x2;
    y->entry[6] = x1 + x3 * x4;
    y->entry[7] = 9 - x1 - x3 * x4;
    y->entry[8] = x4;
    y->entry[9] = 2 - x4;
    y->entry[10] = x5 - 2;
    y->entry[11] = 5 - x5;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 5, 12);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(5, 12, 0, 12, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 5);
    Vector *x0 = vector_alloc(5);

    vector_fill_const(x0, 1.0);

    Vector *lambda0 = vector_alloc(12);
    vector_fill_const(lambda0, 0.0);
    Vector *xstar = vector_alloc(5);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);
    printf(" min f = %lf\n", _fun(xstar));
    return 0;
}
