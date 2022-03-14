#include "matrix.h"
#include "vector.h"
struct constraints
{
    int dim;   // 要约束的空间大小
    int size;  // 约束的数量
    int e;     // 等式数量
    int i;     // 不等式数量 e+i = size; 0,1,2,..,e-1 等式约束,e,e+1,...,e+i-1,大于号约束
    Matrix *A; // 系数矩阵
    Vector *b; // 约束向量
};

typedef struct constraints Constraints;
Constraints *constraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b);
int optimize_qp_linear_constraints(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy);
int optimize_qp_active_set(const Matrix *G, const Vector *c, const Constraints *cons, const Vector *x0, Vector *x_star);
void constrains_subconstrains(const Constraints *con, const Index_set *set, Matrix *A, Vector *b);
void constraints_free(Constraints *con, int recursion);
