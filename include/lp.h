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

LinearConstraints *constraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b);
void constrains_subconstrains(const LinearConstraints *con, const Index_set *set, Matrix *A, Vector *b);
void constraints_free(LinearConstraints *con, int recursion);
void *constraints_verification(const LinearConstraints *con, const Vector *x, Index_set *set);

void optimize_lp_standard_type(const Vector *c, const Vector *b, const Matrix *A, const Vector *x0, Vector *xstar);