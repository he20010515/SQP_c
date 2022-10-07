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

    return 0.2 * (x1 * x1 + x2 * x2 + x3 * x3) + 310 * x1 + 305 * x2 + 300 * x3 - 1000;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    double x3 = x->entry[2];
    y->entry[0] = x1 + x2 + x3 - 240;
    y->entry[1] = x1 - 60;
    y->entry[2] = x1 + x2 - 140;
    y->entry[3] = x2;
    y->entry[4] = x3;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    // proble
    //  min f(x,y,z) = (x-1)^4+(y-2)*z +(z-3)^2
    //  s.t. x^2 + y^2 + z^2 -4 = 0;
    //       x>-0,y>=0,z>=0
    int n = 2;
    int m = 4;
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 5);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 5, 1, 4, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
    Vector *x0 = vector_alloc(3);
    x0->entry[0] = 0;
    x0->entry[1] = 0;
    x0->entry[2] = 0;

    Vector *lambda0 = vector_alloc(5);
    vector_fill_const(lambda0, 0.0);

    Vector *xstar = vector_alloc(3);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
