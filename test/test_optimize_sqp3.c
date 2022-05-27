/*
 * @Author: HeYuwei
 * @Date: 2022-04-26 18:44:11
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-29 19:24:57
 * @FilePath: \SQP_c\test\test_optimize_sqp3.c
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
    double zz = x->entry[2];
    double temp = (xx - 1.0) * (xx - 1.0) * (xx - 1.0) * (xx - 1.0) + (yy - 2.0) * (yy - 2.0) * (yy - 2.0) * (yy - 2.0) + (zz - 3.0) * (zz - 3.0) * (zz - 3.0) * (zz - 3.0);
    return temp;
}

void _c(const Vector *x, Vector *y)
{
    double xx = x->entry[0];
    double yy = x->entry[1];
    double zz = x->entry[2];
    y->entry[0] = xx + yy + zz - 2.0;
    y->entry[1] = -(xx * xx + yy * yy + zz * zz) - 9.0;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    // proble
    //  min f(x,y,z) = (x-1)^4+(y-2)*z +(z-3)^2
    //  s.t. x^2 + y^2 + z^2 -4 = 0;
    //       x>-0,y>=0,z>=0
    int n = 3;
    int m = 2;
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 2);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 2, 1, 1, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
    Vector *x0 = vector_alloc(3);

    x0->entry[0] = 1;
    x0->entry[1] = 2;
    x0->entry[2] = 1;
    Vector *lambda0 = vector_alloc(2);
    lambda0->entry[0] = 0;
    lambda0->entry[1] = 0;

    Vector *xstar = vector_alloc(3);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);

    return 0;
}
