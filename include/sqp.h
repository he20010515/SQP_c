#include "function.h"
#include "matrix.h"
#include "vector.h"
#include "qp.h"

struct nonlinearconstraints
{
    int dim;  //问题维度
    int size; //约束数量
    int e;    //等式约束数量
    int i;    //不等式约束数量
    NdVectorfunction *c;
};
typedef struct nonlinearconstraints Nonlinearconstraints;
void optimize_sqp(const NdsclaFunction *fun,
                  const Nonlinearconstraints *con,
                  const Vector *x0,
                  const Vector *lambda0,
                  Vector *xstar);
Nonlinearconstraints *nonlinearconstraints_alloc(int dim, int size, int e, int i, NdVectorfunction *c);
void nonlinearconstraints_free(Nonlinearconstraints *con);