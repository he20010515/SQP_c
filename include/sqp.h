/*
 * @Author: HeYuwei
 * @Date: 2022-03-25 10:15:12
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-29 17:29:04
 * @FilePath: \SQP_c\include\sqp.h
 * @Description: sqp主要头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "function.h"
#include "matrix.h"
#include "vector.h"
#include "qp.h"

#define SQP_RANDOM_INIT 10

struct nonlinearconstraints
{
    int dim;  // 问题维度
    int size; // 约束数量
    int e;    // 等式约束数量
    int i;    // 不等式约束数量
    NdVectorfunction *c;
};
typedef struct nonlinearconstraints Nonlinearconstraints;
/**
 * @description: sqp算法
 * @param {NdsclaFunction} *fun 目标函数
 * @param {Nonlinearconstraints} *con 非线性约束对象
 * @param {Vector} *x0  起始点
 * @param {Vector} *lambda0 起始点 lambda
 * @param {Vector} *xstar 目标值点
 * @return {*}
 */
void optimize_sqp(const NdsclaFunction *fun, const Nonlinearconstraints *con, const Vector *x0, const Vector *lambda0, Vector *xstar);

/**
 * @description: 申请非线性约束对象的存储空间
 * @param {int} dim 约束维度
 * @param {int} size 约束数量
 * @param {int} e 等式约束数量
 * @param {int} i 不等式约束数量
 * @param {NdVectorfunction} *c n维向量值函数|约束函数
 * @return {*} 指向一个非线性约束对象的指针
 */
Nonlinearconstraints *nonlinearconstraints_alloc(int dim, int size, int e, int i, NdVectorfunction *c);

/**
 * @description: 释放非线性约束对象的存储空间
 * @param {Nonlinearconstraints} *con
 * @return {*}
 */
void nonlinearconstraints_free(Nonlinearconstraints *con);

void optimize_sqpm(const NdsclaFunction *Objective_Fun,
                   const Nonlinearconstraints *con,
                   const Vector *x_init,
                   const Vector *mu_init,
                   const Vector *lam_init,
                   Vector *xstar);
