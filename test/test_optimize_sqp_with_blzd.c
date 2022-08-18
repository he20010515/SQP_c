/*
 * @Author: HeYuwei
 * @Date: 2022-08-11 18:25:25
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-18 12:09:13
 * @FilePath: \SQP_c\test\test_optimize_sqp_with_blzd.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "sqp.h"
#include "blzd_test.h"

#include "windows.h"
typedef double(__stdcall *Fun)(double[]);
Fun fun = NULL;

double _fun(Vector *x)
{
    return fun(x->entry);
}

void _c(const Vector *x, Vector *y)
{
    y->entry[0] = x->entry[0];
    y->entry[1] = x->entry[1];
    y->entry[2] = x->entry[2];
    y->entry[3] = x->entry[3];
    y->entry[4] = x->entry[4];
    y->entry[5] = x->entry[5];
    y->entry[6] = x->entry[6];
    return;
}

int main(int argc, char const *argv[])
{
    sqp_init();

    HMODULE hlib = LoadLibraryA("D:\\Workspace\\HIT\\SQP\\BLZD\\x64\\Release\\BLZD.dll");
    if (hlib == NULL)
    {
        printf("load dll faild");
        return 0;
    }
    fun = (Fun)GetProcAddress(hlib, "_target_function");
    if (fun == NULL)
    {
        printf("getProcAddress fail, error");
    }

    double u0[7] = {0};
    u0[0] = 3.8931935296046505;   //方位角
    u0[1] = 0.0267779154895057;   //极值攻角
    u0[2] = 0.0129239104473724;   //二子级俯仰角初值
    u0[3] = -0.2557582128351035;  //二子级俯仰角终端值
    u0[4] = -0.1745329251994329;  //二子级偏航角初值
    u0[5] = 0.0630259686398771;   //二子级偏航角终端值
    u0[6] = 391.3978486108092056; //二子级末端时间
    double temp;
    Vector *vu0 = vector_alloc(7);
    Vector *lambda0 = vector_alloc(7);
    vector_fill_const(lambda0, 0.0);

    vu0->entry = u0;
    temp = fun(u0);
    printf("targetfunction output:%lf", temp);

    NdVectorfunction *c = ndVectorfunction_alloc(_c, 7, 7);
    Nonlinearconstraints *con = nonlinearconstraints_alloc(7, 7, 0, 7, c);
    NdsclaFunction *f = ndscla_function_alloc(_fun, 7);

    Vector *x0 = vector_alloc(7);
    optimize_sqp(f, con, vu0, lambda0, x0);

    return 0;
}
