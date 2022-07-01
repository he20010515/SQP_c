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
double _fun(Vector *x)
{
    double xx = x->entry[0];
    double yy = x->entry[1];
    double temp = -3.1415926 * xx * xx * yy;
    return temp;
}

void _c(const Vector *x, Vector *y)
{
    double xx = x->entry[0];
    double yy = x->entry[1];
    double zz = x->entry[2];
    y->entry[0] = 3.1415926 * xx * yy + 3.1415926 * xx * xx - 150;
    y->entry[1] = xx;
    y->entry[2] = yy;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    // proble
    //  min f(x,y,z) = (x-1)^4+(y-2)*z +(z-3)^2
    //  s.t. x^2 + y^2 + z^2 -4 = 0;
    //       x>-0,y>=0,z>=0
    int n = 2;
    int m = 3;
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 2, 3);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(2, 3, 1, 2, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 2);
    Vector *x0 = vector_alloc(2);
    x0->entry[0] = 1;
    x0->entry[1] = 2;

    Vector *lambda0 = vector_alloc(3);
    lambda0->entry[0] = 0;
    lambda0->entry[1] = 0;
    lambda0->entry[2] = 0;

    Vector *xstar = vector_alloc(2);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
