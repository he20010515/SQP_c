/*
 * @Author: HeYuwei
 * @Date: 2022-04-22 15:47:42
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-18 19:50:17
 * @FilePath: \SQP_c\test\test_optimize_sqp.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "sqp.h"
#include "math.h"
#include "elog.h"
double _fun(Vector *x) //目标函数定义
{
    double xx = x->entry[0];
    double yy = x->entry[1];
    double zz = x->entry[2];
    double temp = xx * xx - yy * yy - 5 * zz * zz + sin(xx) + cos(yy) - xx * yy * zz - expl(sin(xx) - cos(xx));
    return temp;
}

void _c(const Vector *x, Vector *y) //约束函数定义
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
        // 申请一个NdVectorfunction,为什么要这样做呢? 主要是在结构体中存储关于函数的输入维度以及输出维度,比如本例子中约束函数输入是3维,输出是4维
        NdVectorfunction *c = ndVectorfunction_alloc(_c, 3, 4); 
        // 将申请得到NdVectorfunction与约束的其他信息比如等式约束的数量,不等式约束的数量包装在一个Nonlinearconstraints结构中,可以在头文件中查看函数文档.
        Nonlinearconstraints *con = nonlinearconstraints_alloc(3, 4, 1, 3, c); 
        // 同上,将函数指针包装一下,存储函数的输入维度
        NdsclaFunction *f = ndscla_function_alloc(_fun, 3);
        // 优化起点
        Vector *x0 = vector_alloc(3);
        x0->entry[0] = 3;
        x0->entry[1] = 1;
        x0->entry[2] = 8;
        // 优化起点的lambda乘子
        Vector *lambda0 = vector_alloc(4);
        lambda0->entry[0] = 0;
        lambda0->entry[1] = 0;
        lambda0->entry[2] = 0;
        lambda0->entry[3] = 0;
        // 申请一个向量用于存储最优解
        Vector *xstar = vector_alloc(3);
        // 优化
        optimize_sqp(f, con, x0, lambda0, xstar);
        // 将优化结果打印出来
        vector_print(xstar);
        // 释放存储空间
        ndVectorfunction_free(c);
        nonlinearconstraints_free(con);
        ndscla_function_free(f);
        vector_free(x0);
        vector_free(lambda0);
        vector_free(xstar);
    }

    return 0;
}
