/*
 * @Author: HeYuwei
 * @Date: 2022-03-27 19:10:22
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 18:13:06
 * @FilePath: \SQP_c\include\qp.h
 * @Description:二次优化相关操作头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#pragma once

#include "matrix.h"
#include "vector.h"
#include "lp.h"

/**
 * @description: 求解等式约束下的二次规划问题
 *   min \frac12 x^{T}Hx + c^{T}x
 *  s.t. Ax = b
 * Build Matrix and solve
 *  [H A^T][   x  ] = [-c]
 *  [A  0 ][lambda] = [-b]
 * @param {Matrix} *H 二次型矩阵
 * @param {Vector} *c 一次项向量
 * @param {Matrix} *A 等式约束矩阵
 * @param {Vector} *b 等式约束向量
 * @param {Vector} *x_star 极值点
 * @param {double} *maxy   y的最大值
 * @return {*}
 */
int optimize_qp_linear_constraints(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy);

/**
 * @description: active set方法求解等式以及不等式约束下的二次规划问题
 * 有效集法求解一般约束下的二次优化问题
 * see: https://zhuanlan.zhihu.com/p/29525367
 *
 * mainproblem:
 *  min \frac12 x^{T}Hx + c^{T}x
 * s.t. Ax =(>=) b A (when i = 0,1,..,m  is '=' )
 *                    (when i = m+1,m+2,...,m+k is '>=')
 * @param {Matrix} *G 二次型矩阵
 * @param {Vector} *c 一次项向量
 * @param {LinearConstraints} *cons 线性约束对象
 * @param {Vector} *x0  起始点,若为NULL,则通过lp中的初始点算法计算初始点
 * @param {Vector} *x_star 极值点
 * @param {Vector} *lam 极值点处约束条件的lambda值
 * @return {*}
 */
int optimize_qp_active_set(const Matrix *G, const Vector *c, const LinearConstraints *cons, const Vector *x0, Vector *x_star, Vector *lam);
