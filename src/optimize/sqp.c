#include "matrix.h"
#include "vector.h"
#include "qp.h"
#include "function.h"

struct nonlinearconstraints
{
    int dim;  //问题维度
    int size; //约束数量
    int e;    //等式约束数量
    int i;    //不等式约束数量
    NdVectorfunction *c;
};
typedef struct nonlinearconstraints Nonlinearconstraints;
void optimize_sqp(const NdsclaFunction *fun, const Nonlinearconstraints *con, const Vector *x0, Vector *xstar);
{
    if (__check_input(fun, con, x0, xstar))
    {
        terminate("Please check input shape");
    }

    return;
}

int __check_input(const NdsclaFunction *fun, const Nonlinearconstraints *con, const Vector *x0, const Vector *xstar)
{
    int n = fun->inputSize;
    int m = con->size; // 约束的大小
    int t = !(
        n == fun->inputSize AND n == con->dim AND n == x0->size AND \ 
        n == xstar->size AND\ 
        n == con->c->inputdim AND m == con->size AND m == con->c->outputdim);
    return t;
}
