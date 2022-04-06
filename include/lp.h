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

LinearConstraints *linearconstraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b);
void linearconstrains_subconstrains(const LinearConstraints *con, const Index_set *set, Matrix *A, Vector *b);
void linearconstraints_free(LinearConstraints *con, int recursion);
void *linearconstraints_verification(const LinearConstraints *con, const Vector *x, Index_set *set);

double optimize_lp_standard_type(const Vector *c, const Vector *b, const Matrix *A, const int *init_base, Vector *xstar);
void optimize_lp_2stage(const Vector *c, const Vector *b, const Matrix *A, Vector *xstar);
void optimize_simplex_method(const Vector *_c, const Vector *_b, const Matrix *_A, Vector *xstar, int *init_base);