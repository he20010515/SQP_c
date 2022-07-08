/*
 * @Author: HeYuwei
 * @Date: 2022-04-22 15:47:42
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-08 11:40:06
 * @FilePath: \SQP_c\test\test_optimize_sqp.c
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
    double temp = xx * xx - yy * yy - 5 * zz * zz + sin(xx) + cos(yy) - xx * yy * zz - expl(sin(xx) - cos(xx));
    return temp;
}

void _c(const Vector *x, Vector *y)
{
    y->entry[0] = x->entry[1] + x->entry[2] - 9; // y+z-9 ==0;
    y->entry[1] = x->entry[0] - 1;
    y->entry[2] = x->entry[1] - 1;
    y->entry[3] = x->entry[2] - 7;
}

int main(int argc, char const *argv[])
{
    sqp_init();
    int n = 3;
    int m = 4;

    for (int num = 0; num < 3; num++)
    {
        NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 4);
        Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 4, 1, 3, c);
        NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
        Vector *x0 = vector_alloc(3);

        x0->entry[0] = 3;
        x0->entry[1] = 1;
        x0->entry[2] = 8;
        Vector *lambda0 = vector_alloc(4);
        lambda0->entry[0] = 0;
        lambda0->entry[1] = 0;
        lambda0->entry[2] = 0;
        lambda0->entry[3] = 0;

        Vector *xstar = vector_alloc(3);
        optimize_sqp(f, con, x0, lambda0, xstar);
        vector_print(xstar);

        ndVectorfunction_free(c);
        nonlinearconstraints_free(con);
        ndscla_function_free(f);
        vector_free(x0);
        vector_free(lambda0);
        vector_free(xstar);
    }

    return 0;
}
