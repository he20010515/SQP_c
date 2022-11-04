/*
 * @Author: HeYuwei
 * @Date: 2022-04-26 18:44:11
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-04 22:59:14
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

    return 3 * x1 * x2 + 5 * x3 + 2 * x2 * x4;
}

void _c(const Vector *x, Vector *y)
{
    double x1 = x->entry[0];
    double x2 = x->entry[1];
    double x3 = x->entry[2];
    double x4 = x->entry[3];

    y->entry[0] = x1;
    y->entry[1] = 5 - x1;
    y->entry[2] = -x2 * x3 + 10;
    y->entry[3] = x2 - 3;
    y->entry[4] = x4 - 6;
    y->entry[5] = 15 - x1 - x3 - x4;
    y->entry[6] = x1 * x2 * x3 - 9;
    y->entry[7] = x3;
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
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 4, 8);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(4, 8, 0, 8, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 4);
    Vector *x0 = vector_alloc(4);
    x0->entry[0] = 4;
    x0->entry[1] = 4;
    x0->entry[2] = 4;
    x0->entry[3] = 4;

    Vector *lambda0 = vector_alloc(8);
    vector_fill_const(lambda0, 0.0);
    Vector *xstar = vector_alloc(4);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
