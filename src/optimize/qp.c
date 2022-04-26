#include "matrix.h"
#include "vector.h"
#include "function.h"
#include "util.h"
#include "linear_equations.h"
#include "index_set.h"
#include "qp.h"
#include "math.h"
#include "elog.h"
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define LOG_TAG "qp"
#define MAX_ITER 200

void __qp_compute_subproblem(const Index_set *W_k, const LinearConstraints *cons, const Matrix *G, const Vector *c, Vector *xk, Vector *p, double *y)
{

    int sizeofw = index_set_size(W_k);
    // if (sizeofw > cons->dim)
    //     terminate("ERROR size of wk too big"); //! 如果指标集的大小比输入维度还大的话,这种问题我们目前还没法处理,报错
    //子问题目标函数的一次项
    Vector *Gxk = vector_alloc(cons->dim);
    matrix_mutiply_vector(G, xk, Gxk);
    Vector *Gxk_c = vector_alloc(cons->dim);
    vector_add_vector(Gxk, c, Gxk_c);
    //子问题约束矩阵
    Matrix *sub_A = matrix_alloc(index_set_size(W_k), cons->dim); // 约束个行,dim个列
    Vector *sub_b = vector_alloc(index_set_size(W_k));
    linearconstrains_subconstrains(cons, W_k, sub_A, sub_b);
    vector_fill_const(sub_b, 0);
    //计算子问题得到p
    optimize_qp_linear_constraints(G, Gxk_c, sub_A, sub_b, p, y);
    // matrix_print(G);
    // vector_print(Gxk_c);
    // matrix_print(sub_A);
    // vector_print(sub_b);

    vector_free(Gxk);
    vector_free(Gxk_c);
    matrix_free(sub_A);
    vector_free(sub_b);
    return;
}
Vector *__qp_compute_lambda(const Index_set *W_k, const LinearConstraints *cons, const Matrix *G, const Vector *c, Vector *xk, Index_set *index_set_I, Vector *lambda)
{
    //* 计算拉格朗日系数 lambda i
    // sum_{i\in W}{a_i \lambdai = g = Gx +c}
    // 子线性方程组:
    int m = cons->size;
    Vector *Gxk = vector_alloc(cons->dim);
    matrix_mutiply_vector(G, xk, Gxk);
    Vector *b = vector_alloc(cons->dim);
    vector_add_vector(Gxk, c, b);
    Matrix *R_T = matrix_alloc(index_set_size(W_k), cons->dim);
    Matrix *R = matrix_alloc(cons->dim, index_set_size(W_k));
    matrix_submatrix_by_rowindex_set(cons->A, W_k, R_T);
    matrix_transpose(R_T, R);
    Vector *sub_lambda = vector_alloc(index_set_size(W_k)); // wk大小的lambda

    //求解超定方程组 // min ||Rx-b||
    // x = (R^T *R)^{-1} *RT *y;
    Matrix *RTR = matrix_multiply(R_T, R);
    Matrix *RTR_1 = matrix_alloc(RTR->row_size, RTR->col_size);
    matrix_inverse(RTR, RTR_1);
    Matrix *RTR_1RT = matrix_multiply(RTR_1, R_T);
    matrix_mutiply_vector(RTR_1RT, b, sub_lambda);
    matrix_free(RTR);
    matrix_free(RTR_1);
    matrix_free(RTR_1RT);
    //将lambda 放大到和activeset一样大
    int temp = 0;
    for (int i = 0; i < W_k->index_range; i++)
    {
        if (index_set_is_in(W_k, i))
        {
            lambda->entry[i] = sub_lambda->entry[temp];
            temp++;
        }
        else
            lambda->entry[i] = NAN;
    }
    Index_set *W_k_inter_I = index_set_alloc(m);
    index_set_intersection(W_k, index_set_I, W_k_inter_I);
    Vector *subsublambda = vector_alloc(index_set_size(W_k_inter_I)); // 申请一个和Wk 交 I一样大小的向量,随后将满足条件的lambda 放进去
    temp = 0;
    for (int i = 0; i < m; i++)
    {
        if (index_set_is_in(W_k_inter_I, i))
        {
            subsublambda->entry[temp] = lambda->entry[i];
            temp++;
        }
    }
    vector_free(Gxk);
    vector_free(b);
    matrix_free(R_T);
    matrix_free(R);
    vector_free(sub_lambda);
    index_set_free(W_k_inter_I);
    return subsublambda;
}
double __qp_compute_alphak(const Index_set *W_k, const LinearConstraints *cons, const Vector *x_k, const Vector *p, Vector *alphas)
{
    //计算alphak
    double alphak;
    int m = cons->size;
    int j = 0;
    for (int i = 0; i < m; i++)
    {
        Vector *ai = vector_alloc(cons->A->col_size); // 获取ai
        for (int k = 0; k < cons->A->col_size; k++)
            ai->entry[k] = cons->A->matrix_entry[i][k];
        double aipk = vector_inner_product(ai, p);
        if (aipk < 0 AND !(index_set_is_in(W_k, i)))
        {
            double temp = (cons->b->entry[i] - vector_inner_product(ai, x_k)) / aipk;
            alphas->entry[i] = temp;
        }
        else
            alphas->entry[i] = NAN;
        vector_free(ai);
    }
    alphak = vector_min(alphas);
    return alphak;
}

int optimize_qp_linear_constraints(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy)
{
    // see: https://zhuanlan.zhihu.com/p/375762164
    // see: https://blog.csdn.net/huangdianye/article/details/85030210
    // see: https://blog.csdn.net/qq_42247231/article/details/106938003
    // min \frac12 x^{T}Hx + c^{T}x
    // s.t. Ax = b
    // Build Matrix and solve
    // [H A^T][   x  ] = [-c]
    // [A  0 ][lambda] = [-b]

    //* Check input:
    int n;
    n = H->col_size;
    int m = A->row_size;
    if (!(n == H->row_size AND n == c->size AND n == A->col_size AND A->row_size == b->size))
    {
        // H必须是方阵;一次项;约束矩阵维度匹配;约束数量与b匹配;
        terminate("ERROR: optimize_qp_linear_constraints:input argument does NOT fit, please check input");
    }
    //* alloc workspace
    Matrix *HAAT0 = matrix_alloc(m + n, m + n);
    Vector *xlambda = vector_alloc(m + n);
    Vector *cb = vector_alloc(m + n);
    Vector *xlambda0 = vector_alloc(m + n);
    vector_fill_const(xlambda0, 0);
    //* set big matrix and vector
    //* set HAAT0:
    for (int i = 0; i < m + n; i++)
    {
        for (int j = 0; j < m + n; j++)
        {
            if (i < n AND j < n)
            {
                HAAT0->matrix_entry[i][j] = H->matrix_entry[i][j];
            }
            else if (n <= i AND i < m + n AND j < n)
            {
                // case A
                HAAT0->matrix_entry[i][j] = A->matrix_entry[i - n][j];
            }
            else if (i < n AND n <= j AND j < m + n)
            {
                // case A^T
                HAAT0->matrix_entry[i][j] = A->matrix_entry[j - n][i];
            }
            else if (n <= i AND i < m + n AND n <= j AND j < m + n)
            {
                // case 0
                HAAT0->matrix_entry[i][j] = 0.0;
            }
        }
    }
    //* set cb
    for (int i = 0; i < m + n; i++)
    {
        if (i < n)
        {
            cb->entry[i] = -c->entry[i];
        }
        else
        {
            cb->entry[i] = -b->entry[i - n];
        }
    }
    // // * compute linear function #TODO 这里用一般的迭代法不收敛,目前先用LU分解验证
    //* 现在采用高斯消去法求解线性方程组
    // linear_equation_sor(HAAT0, cb, 0.0001, xlambda0, xlambda, 0.5);
    linear_equation_gaussian_elimination(HAAT0, cb, xlambda);
    // Matrix *inv = matrix_alloc(2 * n, 2 * n);
    // matrix_inverse(HAAT0, inv);
    // matrix_print(inv);
    // matrix_print(HAAT0);
    // matrix_mutiply_vector(HAAT0, cb, xlambda);
    // matrix_free(inv);

    //* split solution
    for (size_t i = 0; i < n; i++)
    {
        x_star->entry[i] = xlambda->entry[i];
    }

    //* free workspace
    matrix_free(HAAT0);
    vector_free(xlambda);
    vector_free(cb);
    vector_free(xlambda0);
    return 1;
}

int optimize_qp_active_set(const Matrix *G, const Vector *c, const LinearConstraints *cons, const Vector *x0, Vector *x_star, Vector *lam)
{
    // 有效集法求解一般约束下的二次优化问题
    // see: https://zhuanlan.zhihu.com/p/29525367

    // mainproblem:
    // min \frac12 x^{T}Hx + c^{T}x
    // s.t. Ax =(>=) b A (when i = 0,1,..,m  is '=' )
    //                   (when i = m+1,m+2,...,m+k is '<=')

    // alloc workspace
    Vector *start_point = NULL;
    int x0_give = true;
    if (x0 == NULL)
    {
        // use simplex method to find a start fesable solution
        start_point = vector_alloc(cons->dim);
        optimize_get_start_feasable_point(cons, (Vector *)start_point, 100, 1e-8, FALSE);
        x0_give = FALSE;
    }
    else
    {
        start_point = (Vector *)x0;
    }

    int m = cons->size;
    Index_set *W_k = index_set_alloc(m); // 工作集
    Index_set *index_set_I = index_set_alloc(m);
    for (int i = cons->e; i < m; i++)
        index_set_append(index_set_I, i);

    Vector *x_k = vector_alloc(cons->dim);
    Vector *x_k_1 = vector_alloc(cons->dim);
    double *y = (double *)malloc(sizeof(double));
    int k = 0;
    //* first compute 首次计算
    // 计算一个初始点以及其对应的工作集
    vector_copy(start_point, x_k);
    linearconstraints_verification(cons, x_k, W_k); // 获取工作集

    // begin compute
    while (1)
    {
        //*计算子问题得到pk
        for (int i = 0; i < cons->e; i++) //确保添加进去等式约束
            if (index_set_is_in(W_k, i))
                continue;
            else
                index_set_append(W_k, i);

        log_i("========iter k = %d=========", k);
        log_i("xk =");
        vector_log(x_k);

        if (vector_have_na(x_k))
        {
            terminate("have nan find");
        }

        Vector *p = vector_alloc(cons->dim);
        __qp_compute_subproblem(W_k, cons, G, c, x_k, p, y);
        if (double_equal(vector_2norm(p), 0.0)) // if p_k = 0
        {
            //*计算lambda
            Vector *lambda = vector_alloc(cons->size);
            Vector *subsublambda = __qp_compute_lambda(W_k, cons, G, c, x_k, index_set_I, lambda);
            log_i("lambda:");
            // vector_print(lambda);
            if (vector_any_bigger_equal_than_const(subsublambda, 0)) //若lambda i >0 (激活不等式约束集) (\any i \in Wk)
            {
                log_i("case: iter done");
                vector_copy(x_k, x_star);
                // vector_print(x_star);
                vector_copy(lambda, lam);
                vector_free(lambda);
                vector_free(subsublambda);
                vector_free(p);
                break;
            }
            else
            {
                log_i("case: remove con");
                int j = vector_argmin(lambda);
                log_i("remove cons %d", j);
                index_set_remove(W_k, j);
                vector_copy(x_k, x_k_1);
                vector_free(lambda);
                vector_free(subsublambda);
            }
        }
        else
        {
            //*计算alphak
            //更新xk
            Vector *alphas = vector_alloc(cons->size);
            double alphak = 0.0;
            alphak = __qp_compute_alphak(W_k, cons, x_k, p, alphas);
            if (isnan(alphak))
                alphak = 1;
            alphak = min(1, alphak);
            Vector *alphapk = vector_multiply_const(p, alphak, 1); // copy = 1
            int j = vector_argmin(alphas);
            vector_add_vector(x_k, alphapk, x_k_1);
            log_i("alpha k:%lf", alphak);
            // vector_print(alphas);
            vector_free(alphapk);
            vector_free(alphas);
            if (alphak < 1.) //若不满足某些约束
            {
                log_i("case: update wk,and update wk");
                index_set_append(W_k, j); //将不满足的约束添加进工作集
                log_i("Index set append%d", j);
                log_i("alphak = %f", alphak);
            }
            else
            {
                log_i("case: update xk and keep wk"); //约束集不变
                log_i("alphak = %f", alphak);
            }
        }
        log_i("iter%d done;xk = ", k);
        // vector_print(x_k_1);
        k++;
        if (k >= MAX_ITER)
        {
            log_e("Iteration overflow ,Please check inputs. qp is Break");
            break;
        }
        vector_copy(x_k_1, x_k);
        vector_free(p);
    }

    // free workspace
    index_set_free(W_k);
    index_set_free(index_set_I);
    vector_free(x_k);
    vector_free(x_k_1);
    free(y);
    if (!x0_give)
        vector_free((Vector *)start_point);

    return 0;
}
