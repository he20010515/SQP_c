#include "matrix.h"
#include "vector.h"
#include "qp.h"
#include "function.h"
#include "sqp.h"
#include "elog.h"
#include <math.h>

#define LOG_TAG "SQP"

#define NUMERICAL_DIFF_STEP 0.001

struct __wrapper_info
{
    const NdsclaFunction *fun;
    const Nonlinearconstraints *con;
    Vector *lambda;
};

struct __wrapper_info _wrapper_info;

int __check_input(const NdsclaFunction *fun, const Nonlinearconstraints *con, const Vector *x0, const Vector *xstar);
int __check_inner_loop(const Vector *xk, const Vector *pk, double aita, double miu, double alphak, const Nonlinearconstraints *con, const NdsclaFunction *fun);
double __phi_1(const Vector *x, const double miu, const Nonlinearconstraints *con, const NdsclaFunction *fun);
double __D(const Vector *xk, const Vector *pk, const double miu, const NdsclaFunction *fun, const Nonlinearconstraints *con);
double __lagrange_function(Vector *xk, const Vector *lambda, const NdsclaFunction *fun, const Nonlinearconstraints *con);
double __lagrange_wrapper(Vector *x);
void __BFGS_update(const Matrix *Bk, const Vector *lambdak_1, const Vector *xk, const Vector *xk_1, Matrix *Bk_1, NdsclaFunction *lagrange);
// TODO List
//*1. BFGS 修正 // Hession矩阵计算
// 2. lambda 到底应该是多少
// 3. 子问题初始可行解

void optimize_sqp(const NdsclaFunction *fun,
                  const Nonlinearconstraints *con,
                  const Vector *x0,
                  const Vector *lambda0,
                  Vector *xstar)
{
    int iternum = 0;
    if (__check_input(fun, con, x0, xstar))
    {
        terminate("Please check input shape");
    }
    // basic problem info:
    int n = x0->size;
    int m = con->c->outputdim; // 约束数量

    // alloc workspace
    double f0;
    const double rho = 0.5;
    const double aita = 0.25;
    const double tao = 0.5;

    Vector *xk = vector_alloc(n);
    Vector *xk_1 = vector_alloc(n);

    Vector *p = vector_alloc(n);

    Vector *gradf0 = vector_alloc(n);
    Vector *gradfk = vector_alloc(n);
    Vector *gradfk_1 = vector_alloc(n);

    Matrix *A0 = matrix_alloc(m, n); // TODO ensure size of jacobi c;
    Matrix *Ak = matrix_alloc(m, n);
    Matrix *Ak_1 = matrix_alloc(m, n);

    Vector *c0 = vector_alloc(m);
    Vector *ck = vector_alloc(m);
    Vector *ck_1 = vector_alloc(m);

    // Matrix *HxxL0 = matrix_alloc(n, n);
    // Matrix *HxxLk = matrix_alloc(n, n);
    // Matrix *HxxLk_1 = matrix_alloc(n, n);
    Matrix *B0 = matrix_alloc(n, n);
    Matrix *Bk = matrix_alloc(n, n);
    Matrix *Bk_1 = matrix_alloc(n, n);

    Vector *lambdak = vector_alloc(m); //*传给拉格朗日函数的lambda
    Vector *lambdak_1 = vector_alloc(m);
    Vector *plambda = vector_alloc(m); //* lambda 步长
    Vector *lambdahat = vector_alloc(m);

    // init varlable
    vector_copy(x0, xk);
    vector_copy(x0, xk_1);

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
    _wrapper_info.lambda = lambdak;
    vector_copy(lambda0, lambdak);
    NdsclaFunction *lagrange_function = ndscla_function_alloc(__lagrange_wrapper, n);

    matrix_fill_const(B0, 0.0);
    for (int i = 0; i < n; i++)
        B0->matrix_entry[i][i] = 1.;
    matrix_copy(B0, Bk);
    matrix_copy(B0, Bk_1);

    log_i("x0");
    vector_log(x0);
    log_i("Hession matrix of lagrange function");
    matrix_log(B0);
    log_i("grad of target function in x0");
    vector_log(gradf0);
    log_i("jacobian matrix of c in x0");
    matrix_log(A0);

    int k = 0;
    //  mainloop
    while (1)
    {
        // 计算子问题
        Vector *_ck = vector_multiply_const(ck, -1., 1);
        LinearConstraints *subcon = linearconstraints_alloc(n, m, con->e, con->i, Ak, _ck);
        log_i("=========iter k = %d =========", k);
        log_i("Xk = ");
        vector_log(xk);
        log_i("subproblem :");
        log_i("B:");
        matrix_log(Bk);
        log_i("gradfk");
        vector_log(gradfk);
        log_i("subcon:");
        matrix_print(subcon->A);
        vector_print(subcon->b);
        log_i("subproble start point");
        vector_log(xk_1);
        optimize_qp_active_set(Bk, gradfk, subcon, xk_1, p, lambdahat); // 用上一步的结束值当做这一阶段的初值
        log_i("subproblem ans P:");
        vector_log(p);
        log_i("subproblem lambdahat");
        vector_fillna(lambdahat);
        vector_log(lambdahat);

        vector_free(_ck);

        // plambda = lambdahat-lambdak
        Vector *_lambdak = vector_multiply_const(lambdak, -1., 1);
        vector_add_vector(lambdahat, _lambdak, plambda);
        vector_free(_lambdak);

        // choose miuk
        Vector *temp = vector_alloc(n);
        vector_mutiply_matrix(p, Bk, temp);
        double miu = (vector_inner_product(gradfk, p) + 0.5 * vector_inner_product(temp, p)) / ((1 - rho) * vector_1norm(ck));
        vector_free(temp);
        double alphak = 1;
        while (__check_inner_loop(xk, p, aita, miu, alphak, con, fun))
        {
            // update alphak
            alphak = alphak * tao;
        }
        // update xk,lambdak
        vector_add_vector(xk, p, xk_1);

        log_i("xk_1:");
        vector_log(xk_1);
        vector_add_vector(lambdak, plambda, lambdak_1);
        log_i("lambdak_1");
        vector_log(lambdak_1);

        ndscla_central_grad(fun, NUMERICAL_DIFF_STEP, xk_1, gradfk_1);
        ndVectorfunction_call(con->c, xk_1, ck_1);
        ndVectorfunction_jacobian(con->c, xk_1, NUMERICAL_DIFF_STEP, Ak_1);
        __BFGS_update(Bk, _lambdak, xk, xk_1, Bk_1, lagrange_function);
        log_i("Bk_1");
        matrix_print(Bk_1);

        // swap
        vector_copy(xk_1, xk);
        vector_copy(gradfk_1, gradfk);
        vector_copy(ck_1, ck);
        matrix_copy(Ak_1, Ak);
        matrix_copy(Bk_1, Bk);
        iternum++;
        if (iternum >= 100)
        {
            log_e("iter overflow");
            break;
        }
    }
    // // free workspace

    return;
}

Nonlinearconstraints *nonlinearconstraints_alloc(int dim, int size, int e, int i, NdVectorfunction *c)
{
    Nonlinearconstraints *con = (Nonlinearconstraints *)malloc(sizeof(Nonlinearconstraints));
    con->dim = dim;
    con->size = size;
    con->e = e;
    con->i = i;
    con->c = c;
    return con;
}

void nonlinearconstraints_free(Nonlinearconstraints *con)
{
    free(con);
}

void __BFGS_update(const Matrix *Bk, const Vector *lambdak_1, const Vector *xk, const Vector *xk_1, Matrix *Bk_1, NdsclaFunction *lagrange)
{
    Vector *sk = vector_alloc(xk->size);
    Vector *_xk = vector_multiply_const(xk, -1., 1);
    vector_add_vector(xk_1, _xk, sk); // sk = xK+1-xk;

    Vector *gradxk_1 = vector_alloc(xk->size);
    ndscla_central_grad(lagrange, NUMERICAL_DIFF_STEP, xk, gradxk_1);
    Vector *gradxk = vector_alloc(xk->size);
    ndscla_central_grad(lagrange, NUMERICAL_DIFF_STEP, xk_1, gradxk);
    Vector *yk = vector_alloc(xk->size);
    Vector *_gradxk = vector_multiply_const(gradxk, -1., 1);
    vector_add_vector(gradxk_1, _gradxk, yk);

    // thetak
    double thetak = 0;
    double skyk = vector_inner_product(sk, yk);
    Vector *temp = vector_alloc(Bk->row_size);
    vector_mutiply_matrix(sk, Bk, temp);
    double SBS = vector_inner_product(temp, sk);
    vector_free(temp);
    if (skyk >= SBS)
        thetak = 1;
    else
        thetak = SBS * 0.8 / (SBS - skyk);

    // compute r r= \theta y + (1-theta)Bksk
    Vector *thetay = vector_multiply_const(yk, thetak, 1);
    Vector *Bksk = vector_alloc(Bk->col_size);
    matrix_mutiply_vector(Bk, sk, Bksk);
    Vector *one_min_theta_bksk = vector_multiply_const(Bksk, (1. - thetak), 1);
    Vector *r = vector_alloc(yk->size);
    vector_add_vector(thetay, one_min_theta_bksk, r);

    // compute \frac{Bss^TB}{s^TBs}
    Matrix *ssT = matrix_alloc(sk->size, sk->size);
    vector_mutiply_vectorT(sk, sk, ssT);
    Matrix *BssT = matrix_multiply(Bk, ssT);
    Matrix *BssTB = matrix_multiply(BssT, Bk);
    Matrix *_BssTB_SBS = matrix_alloc(BssT->row_size, BssT->col_size);
    matrix_mutiply_const(BssTB, -1. / SBS, _BssTB_SBS);

    // compute \frac{rr^T}{s^Tr}
    // Matrix *rrT = matrix_alloc();
    Matrix *rrt = matrix_alloc(Bk->col_size, Bk->col_size);
    vector_mutiply_vectorT(r, r, rrt);
    Matrix *rrt_st = matrix_alloc(rrt->col_size, rrt->col_size);
    matrix_mutiply_const(rrt, 1 / vector_inner_product(sk, r), rrt_st);
    Matrix *tempmat = matrix_alloc(Bk->col_size, Bk->row_size);
    //*update Bk
    matrix_add(tempmat, _BssTB_SBS, rrt_st);
    matrix_add(Bk_1, (Matrix *)Bk, (Matrix *)tempmat);

    vector_free(sk);
    vector_free(_xk);
    vector_free(gradxk);
    vector_free(gradxk_1);
    vector_free(_gradxk);
    vector_free(yk);

    vector_free(thetay);
    vector_free(Bksk);
    vector_free(one_min_theta_bksk);

    vector_free(r);
    matrix_free(rrt);
    matrix_free(rrt_st);
    matrix_free(tempmat);

    matrix_free(ssT);
    matrix_free(BssT);
    matrix_free(BssTB);
    matrix_free(_BssTB_SBS);
}

int __check_input(const NdsclaFunction *fun, const Nonlinearconstraints *con, const Vector *x0, const Vector *xstar)
{
    int n = fun->inputSize;
    int m = con->size; // 约束的大小
    int t = !(
        n == fun->inputSize AND n == con->dim AND n == x0->size AND n == xstar->size AND n == con->c->inputdim AND m == con->size AND m == con->c->outputdim);
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

    log_i("abs(l-r)%lf", fabs(l - r));
    return l > r;
}

double __phi_1(const Vector *x, const double miu, const Nonlinearconstraints *con, const NdsclaFunction *fun)
{
    Vector *cx = vector_alloc(con->c->outputdim);
    ndVectorfunction_call(con->c, x, cx);
    Vector *xw = (Vector *)x;
    double fx = ndscla_function_call(fun, xw);
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

double __lagrange_function(Vector *xk, const Vector *lambda, const NdsclaFunction *fun, const Nonlinearconstraints *con)
{
    Vector *cx = vector_alloc(con->c->outputdim);
    ndVectorfunction_call(con->c, xk, cx);
    return ndscla_function_call(fun, xk) - vector_inner_product(lambda, cx);
}

double __lagrange_wrapper(Vector *x)
{
    // 借助一个文件内的结构体_wrapper_info,将拉格朗日函数包装为可以用Ndsclafunciton Hession可以求导的样式
    return __lagrange_function(x, _wrapper_info.lambda, _wrapper_info.fun, _wrapper_info.con);
}