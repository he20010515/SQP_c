#pragma once
#include "vector.h"
#include "matrix.h"

struct simplex
{
    Vector *c;    // 系数向量
    Matrix *A_ub; // 不等式约束矩阵
    Vector *b_ub; // 不等式约束边界
    Matrix *A_eq; // 等式约束矩阵
    Vector *b_eq; // 等式约束条件
    int N_x;      // 变量个数
    int meq;      // 等式约束个数
    int mub;      // 不等式约束个数
    int m;        // 约束总数
    Matrix *T;    // 单纯形表
    int *sign;    // 记录基变量序号
    int F;        // FLAG量 1,继续迭代,0,无界解,停止迭代
};
typedef struct simplex Simplex;
Simplex *simplex_alloc(const Vector *c, const Matrix *A_ub, const Vector *b_ub, const Matrix *A_eq, const Vector *b_eq);
void simplex_inital_value(Simplex *self);
void simplex_free(Simplex *self);
void simplex_solve(Simplex *self);
int simplex_calculate(Simplex *self);
void simplex_change(Simplex *self);
int simplex_main(Simplex *self, Vector *xstar);
