/*
 * @Author: HeYuwei
 * @Date: 2022-04-01 18:30:23
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-18 12:10:17
 * @FilePath: \SQP_c\test\test_optimize_sqp2.c
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
    double temp = (xx - 1) * (xx - 1) * (xx - 1) * (xx - 1) + (yy - 2) * zz + (zz - 3) * (zz - 3);
    return temp;
}

void _c(const Vector *x, Vector *y)
{
    y->entry[0] = x->entry[0] - 4;
    y->entry[1] = x->entry[0];
    y->entry[2] = x->entry[1];
    y->entry[3] = x->entry[2];
}

int main(int argc, char const *argv[])
{
    sqp_init();
    // proble
    //  min f(x,y,z) = (x-1)^4+(y-2)*z +(z-3)^2
    //  s.t. x^2 + y^2 + z^2 -4 = 0;
    //       x>-0,y>=0,z>=0
    int n = 3;
    int m = 4;
    NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 4);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 4, 1, 3, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
    Vector *x0 = vector_alloc(3);

    x0->entry[0] = 0;
    x0->entry[1] = 0;
    x0->entry[2] = sqrtl(2);
    Vector *lambda0 = vector_alloc(4);
    lambda0->entry[0] = 0;
    lambda0->entry[1] = 0;
    lambda0->entry[2] = 0;
    lambda0->entry[3] = 0;

    Vector *xstar = vector_alloc(3);
    optimize_sqp(f, con, x0, lambda0, xstar);
    vector_print(xstar);
    log_i("number of targetfunction call %ld times", f->call_num);

    return 0;
}
