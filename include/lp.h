/*
 * @Author: HeYuwei
 * @Date: 2022-04-03 19:38:47
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 11:12:04
 * @FilePath: \SQP_c\include\lp.h
 * @Description: 线性规划求解头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#pragma once

#include "vector.h"
#include "matrix.h"

struct linearconstraints
{
    int dim;   // 要约束的空间大小
    int size;  // 约束的数量
    int e;     // 等式数量
    int i;     // 不等式数量 e+i = size; 0,1,2,..,e-1 等式约束,e,e+1,...,e+i-1,大于号约束
    Matrix *A; // 系数矩阵
    Vector *b; // 约束向量
};
typedef struct linearconstraints LinearConstraints;

/**
 * @description: 申请线性约束结构体
 * @param {int} dim   问题维度
 * @param {int} size  约束数量
 * @param {int} e     等式约束数量
 * @param {int} i     大于等于约束数量
 * @param {Matrix} *A Ax  = b
 * @param {Vector} *b Ax >= b
 * @return {*}
 */
LinearConstraints *linearconstraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b);

/**
 * @description: 根据指标器Index_set建立自约束  通过给定的指标集和指定的约束构建子等式约束  A矩阵应为与G相等的方阵,b应与其匹配的列向量
 * @param {LinearConstraints} *con 约束
 * @param {Index_set} *set 指标集
 * @param {Matrix} *A
 * @param {Vector} *b
 * @return {*}
 */
void linearconstrains_subconstrains(const LinearConstraints *con, const Index_set *set, Matrix *A, Vector *b);

/**
 * @description: 释放线性约束结构体存储空间
 * @param {LinearConstraints} *con
 * @param {int} recursion 是否递归释放,若为true,则释放其指向的A以及b
 * @return {*}
 */
void linearconstraints_free(LinearConstraints *con, int recursion);

/**
 * @description: 验证当前点x所满足的约束(全部视为等式约束),输出为指标集
 * @param {LinearConstraints} *con 约束
 * @param {Vector} *x 当前点
 * @param {Index_set} *set 指标集
 * @return {*}
 */
void linearconstraints_verification(const LinearConstraints *con, const Vector *x, Index_set *set);

/**
 * @description: 求解线性规划
 *       min c^T x
 *  s.t. Ax == b (i < e    )
 *       Ax >= b (e <=i < m)
 * @param {LinearConstraints} *con 线性约束
 * @param {Vector} *c 约束指标
 * @param {Vector} *x 解
 * @param {int} maxiter 最大迭代次数
 * @param {double} tol 容差
 * @param {int} bland  是否使用bland法则
 * @return {*}
 */
int optimize_lp(const LinearConstraints *con, const Vector *c, Vector *x, int maxiter, double tol, int bland);

/**
 * @description: 求解线性约束con下的一个初始可行点
 * @param {LinearConstraints} *con 线性约束
 * @param {Vector} *x0  可行点
 * @param {int} maxiter 最大迭代次数
 * @param {double} tol  容差
 * @param {int} bland   是否使用bland法则
 * @return {*}
 */
int optimize_get_start_feasable_point(const LinearConstraints *con, Vector *x0, int maxiter, double tol, int bland);