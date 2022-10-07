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

    return 1.5 * x1 * x1 + x2 * x2 + 0.85 * x3 * x3 + 3 * x1 - 8.2 * x2 - 1.95 * x3;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    double x3 = x->entry[2];

    y->entry[0] = x1 + x2 + x3 - 3;
    y->entry[1] = 2 - x1 - x3;
    y->entry[2] = 2 + x1 - 2 * x2;
    y->entry[3] = 3 - x2 - 2 * x3;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 4);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 4, 1, 3, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
    Vector *x0 = vector_alloc(3);
    x0->entry[0] = 0.1;
    x0->entry[1] = 1.1;
    x0->entry[2] = 0.1;

    Vector *lambda0 = vector_alloc(4);
    vector_fill_const(lambda0, 0.0);
    Vector *xstar = vector_alloc(3);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
