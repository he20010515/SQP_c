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

    return 2 * x1 * x1 - 4 * x1 * x2 + 4 * x2 * x2 - 6 * x1 - 3 * x2;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    y->entry[0] = -x1 - x2 + 3;
    y->entry[1] = -4 * x1 - x2 + 9;
    y->entry[2] = x1;
    y->entry[3] = x2;
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
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 2, 4);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(2, 4, 0, 4, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 2);
    Vector *x0 = vector_alloc(2);
    x0->entry[0] = 0;
    x0->entry[1] = 0;

    Vector *lambda0 = vector_alloc(4);
    lambda0->entry[0] = 0;
    lambda0->entry[1] = 0;
    lambda0->entry[2] = 0;

    Vector *xstar = vector_alloc(2);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
