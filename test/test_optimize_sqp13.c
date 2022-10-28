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

    return x1 * x2 * x2;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];

    y->entry[0] = x1 * x1 + x2 * x2 - 2;
    return;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 2, 1);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(2, 1, 1, 0, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 2);
    Vector *x0 = vector_alloc(2);
    Vector *xstar = vector_alloc(2);
    Vector *lambda0 = vector_alloc(1);

    vector_fill_const(x0, 1.0);
    vector_fill_const(lambda0, 0.0);

    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);
    printf(" min f = %lf\n", _fun(xstar));
    return 0;
}
