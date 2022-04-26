/*
 * @Author: HeYuwei
 * @Date: 2022-03-13 09:13:40
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 12:09:51
 * @FilePath: \SQP_c\include\linear_equations.h
 * @Description: 线性方程组求解
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "vector.h"
#include "matrix.h"
#pragma once
/**
 * @description: 雅克比迭代法求解线性方程组  Ax=b
 * @param {Matrix} *mat A矩阵
 * @param {Vector} *b   b向量
 * @param {double} epsilon 容差
 * @param {Vector} *x0  起始点
 * @param {Vector} *x   解
 * @return {*}
 */
void linear_equation_jacobi(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x);
/**
 * @description:  高斯赛德尔迭代求解线性方程组 Ax=b
 * @param {Matrix} *mat A矩阵
 * @param {Vector} *b   b向量
 * @param {double} epsilon 容差
 * @param {Vector} *x0 起始点
 * @param {Vector} *x  解
 * @return {*}
 */
void linear_equation_gauss_sidel(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x);

/**
 * @description: 超松弛迭代求解线性方程组 Ax = b
 * @param {Matrix} *mat A矩阵
 * @param {Vector} *b   b向量
 * @param {double} epsilon 容差
 * @param {Vector} *x0  起始点
 * @param {Vector} *x   解
 * @param {double} omega 超松弛迭代系数
 * @return {*}
 */
void linear_equation_sor(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x, const double omega);

/**
 * @description: 高斯消去法求解线性方程组 Ax = b
 * @param {Matrix} *mat A 矩阵
 * @param {Vector} *b   b 向量
 * @param {Vector} *x   x 解
 * @return {*}
 */
void linear_equation_gaussian_elimination(const Matrix *mat, const Vector *b, Vector *x);