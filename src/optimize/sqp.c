#include "matrix.h"
#include "vector.h"
#include "qp.h"
#include "function.h"

#define NUMERICAL_DIFF_STEP 0.001

struct nonlinearconstraints
{
    int dim;  //问题维度
    int size; //约束数量
    int e;    //等式约束数量
    int i;    //不等式约束数量
    NdVectorfunction *c;
};
typedef struct nonlinearconstraints Nonlinearconstraints;

struct __wrapper_info
{
    NdsclaFunction *fun;
    Nonlinearconstraints *con;
    Vector *lambda;
};

struct __wrapper_info _wrapper_info;
void optimize_sqp(const NdsclaFunction *fun,
                  const Nonlinearconstraints *con,
                  const Vector *x0,
                  const Vector *lambda0,
                  Vector *xstar)
{
    if (__check_input(fun, con, x0, xstar))
    {
        terminate("Please check input shape");
    }
    // basic problem info:
    int n = x0->size;
    int m = con->c->outputdim; // 约束数量

    // alloc workspace
    double f0;
    Vector *gradf0 = vector_alloc(n);
    Vector *gradfk = vector_alloc(n);
    Vector *gradfk_1 = vector_alloc(n);

    Matrix *A0 = matrix_alloc(m, n); // TODO ensure size of jacobi c;
    Matrix *Ak = matrix_alloc(m, n);
    Matrix *Ak_1 = matrix_alloc(m, n);

    Vector *c0 = vector_alloc(m);
    Vector *ck = vector_alloc(m);
    Vector *ck_1 = vector_alloc(m);

    Matrix *HxxL0 = matrix_alloc(n, n);
    Matrix *HxxLk = matrix_alloc(n, n);
    Matrix *HxxLk_1 = matrix_alloc(n, n);

    // init varlable
    ndscla_central_grad(fun, NUMERICAL_DIFF_STEP, x0, gradf0);
    vector_copy(gradf0, gradfk);
    vector_copy(gradf0, gradfk_1);

    ndVectorfunction_jacobian(con->c, x0, NUMERICAL_DIFF_STEP, A0);
    matrix_copy(A0, Ak);
    matrix_copy(A0, Ak_1);

    ndVectorfunction_call(con->c, x0, c0);
    vector_copy(c0, ck);
    vector_copy(c0, ck_1);

    _wrapper_info.con = con;
    _wrapper_info.fun = fun;
    _wrapper_info.lambda = lambda0;
    ndscla_central_hession(__lagrange_wrapper, NUMERICAL_DIFF_STEP, x0, HxxL0);
    matrix_copy(HxxL0, HxxLk);
    matrix_copy(HxxL0, HxxLk_1);

    

    //  mainloop
    while (1)
    {
        //计算子问题
        // plambda = lambdahat-lambdak
        // chose miuk
        // alpha = 1
        while (1)
        {
            // update alphak
        }
        // update xk,lambdak
        //构建子问题
        //  如果采用拟牛顿近似 更新bk
    }
    // free workspace

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

int __check_inner_loop(const Vector *xk, const Vector *pk, double aita, double miu, double alphak, const Nonlinearconstraints *con, const NdsclaFunction *fun)
{
    Vector *xk_add_alphapk = vector_alloc(xk->size);
    Vector *alphapk = vector_multiply_const(pk, alphak, 1);
    vector_add_vector(alphapk, xk, xk_add_alphapk);
    double l = __phi_1(xk_add_alphapk, miu, con, fun);
    double r = __phi_1(xk, miu, con, fun) + aita * __D(xk, pk, miu, fun, con);
    vector_free(alphapk);
    vector_free(xk_add_alphapk);
    return l > r;
}

double __phi_1(const Vector *x, const double miu, const Nonlinearconstraints *con, const NdsclaFunction *fun)
{
    Vector *cx = vector_alloc(con->c->outputdim);
    ndVectorfunction_call(con->c, x, cx);
    double fx = ndscla_function_call(fun, x);
    return fx + miu * vector_1norm(cx);
}

double __D(const Vector *xk, const Vector *pk, const double miu, const NdsclaFunction *fun, const Nonlinearconstraints *con)
{
    Vector *gradfk = vector_alloc(fun->inputSize);
    ndscla_central_grad(fun, NUMERICAL_DIFF_STEP, xk, gradfk);
    Vector *ck = vector_alloc(con->c->outputdim);
    ndVectorfunction_call(con->c, xk, ck);
    double temp = vector_inner_product(gradfk, pk) - miu * vector_1norm(ck);
    vector_free(gradfk);
    vector_free(ck);
    return temp;
}

double __lagrange_function(const Vector *xk, const Vector *lambda, const NdsclaFunction *fun, const Nonlinearconstraints *con)
{
    Vector *cx = vector_alloc(con->c->outputdim);
    ndVectorfunction_call(con->c, xk, cx);
    return ndscla_function_call(fun, xk) - vector_inner_product(lambda, cx);
}

double __lagrange_wrapper(const Vector *x)
{
    // 借助一个文件内的结构体_wrapper_info,将拉格朗日函数包装为可以用Ndsclafunciton Hession可以求导的样式
    return __lagrange_function(x, _wrapper_info.lambda, _wrapper_info.fun, _wrapper_info.con);
}