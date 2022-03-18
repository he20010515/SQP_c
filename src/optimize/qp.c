#include "matrix.h"
#include "vector.h"
#include "function.h"
#include "util.h"
#include "linear_equations.h"
#include "index_set.h"
#include "qp.h"

Constraints *constraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b)
{
    Constraints *con = (Constraints *)malloc(sizeof(Constraints));
    if (!(size == e + i AND dim == A->col_size AND size == A->row_size))
    {
        terminate("ERROR :constraints_alloc: numbers of constraints doesn't fit");
    }
    con->dim = dim;
    con->size = size;
    con->e = e;
    con->i = i;
    con->A = A;
    con->b = b;
}

void *constraints_verification(const Constraints *con, const Vector *x, Index_set *set)
{
    // 验证解x是否满足约束,将满足的约束编号添加到约束集合mat上
    if (!(con->dim == x->size AND set->index_range == con->A->row_size))
    {
        terminate("ERROR : constraints_verification: dim of vector to verificate doesn't fit");
    }
    Vector *b_ = vector_alloc(con->b->size);
    matrix_mutiply_vector(con->A, x, b_);
    for (int i = 0; i < b_->size; i++)
    {
        if (i < con->e)
        {
            if (double_equal(b_->entry[i], con->b->entry[i]))
            {
                index_set_append(set, i); // 满足等式约束条件则添加
            }
        }
        else
        {
            if (double_equal(b_->entry[i], con->b->entry[i]))
            {
                index_set_append(set, i); // 满足不等式约束条件则添加
            }
        }
    }
    vector_free(b_);
}

void constraints_free(Constraints *con, int recursion)
{
    if (recursion)
    {
        matrix_free(con->A);
        vector_free(con->b);
        free(con);
    }
    else
    {
        free(con);
    }
}

void constrains_subconstrains(const Constraints *con, const Index_set *set, Matrix *A, Vector *b)
{
    // 通过给定的指标集和指定的约束构建子等式约束
    // A矩阵应为与G相等的方阵,b应与其匹配的列向量
    if (!(con->size == set->index_range AND A->row_size == A->col_size AND A->col_size == con->dim) AND b->size == con->dim)
    {
        terminate("ERROR constrains_subconstrains size not fit");
    }
    int sub_i = 0;
    matrix_fill_const(A, 0.0);
    vector_fill_const(b, 0.0);
    for (int i = 0; i < con->size; i++)
    {
        if (index_set_is_in(set, i))
        {
            // 如果在指标集中
            for (int j = 0; j < A->col_size; j++) // 复制矩阵
            {
                A->matrix_entry[sub_i][j] = con->A->matrix_entry[i][j];
            }
            b->entry[sub_i] = con->b->entry[i];
            sub_i += 1;
        }
        else
            continue;
    }
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
    if (!(n == H->row_size AND n == c->size AND n == A->row_size AND n == A->col_size AND n == b->size))
    {
        terminate("ERROR: optimize_qp_linear_constraints:input argument does NOT fit, please check input");
    }
    //* alloc workspace
    Matrix *HAAT0 = matrix_alloc(2 * n, 2 * n);
    Vector *xlambda = vector_alloc(2 * n);
    Vector *cb = vector_alloc(2 * n);
    Vector *xlambda0 = vector_alloc(2 * n);
    vector_fill_const(xlambda0, 0);
    //* set big matrix and vector
    //* set HAAT0:
    for (int i = 0; i < 2 * n; i++)
    {
        for (int j = 0; j < 2 * n; j++)
        {
            if (i < n AND j < n)
            {
                HAAT0->matrix_entry[i][j] = H->matrix_entry[i][j];
            }
            else if (n <= i AND i < 2 * n AND j < n)
            {
                // case A
                HAAT0->matrix_entry[i][j] = A->matrix_entry[i - n][j];
            }
            else if (i < n AND n <= j AND j < 2 * n)
            {
                // case A^T
                HAAT0->matrix_entry[i][j] = A->matrix_entry[j - n][i];
            }
            else if (n <= i AND i < 2 * n AND n <= j AND j < 2 * n)
            {
                // case 0
                HAAT0->matrix_entry[i][j] = 0.0;
            }
        }
    }
    //* set cb
    for (int i = 0; i < 2 * n; i++)
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

int optimize_qp_active_set(const Matrix *G, const Vector *c, const Constraints *cons, const Vector *x0, Vector *x_star)
{
    // 有效集法求解一般约束下的二次优化问题
    // see: https://zhuanlan.zhihu.com/p/29525367

    // mainproblem:
    // min \frac12 x^{T}Hx + c^{T}x
    // s.t. Ax =(<=) b A (when i = 0,1,..,m  is '=' )
    //                   (when i = m+1,m+2,...,m+k is '<=')

    // alloc workspace
    int m = cons->size;
    Index_set *W_k = index_set_alloc(m); // 工作集
    Index_set *W_k_1 = index_set_alloc(m);
    Vector *x_k = vector_alloc(cons->dim);
    Vector *x_k_1 = vector_alloc(cons->dim);
    double *y = (double *)malloc(sizeof(double));
    int k = 0;
    //* first compute 首次计算
    // 计算一个初始点以及其对应的工作集
    vector_copy(x0, x_k);
    constraints_verification(cons, x_k, W_k); // 获取工作集
    // begin compute
    while (1)
    {
        //*计算子问题得到pk
        int sizeofw = index_set_size(W_k); // TODO 将等式约束加入指标集  DONE
        for (int i = 0; i < cons->e; i++)  //确保添加进去等式约束
            index_set_append(W_k, i);
        printf("work set of x0:\n");
        index_set_print(W_k);
        if (sizeofw > cons->dim)
            terminate("ERROR size of wk too big"); //! 如果指标集的大小比输入维度还大的话,这种问题我们目前还没法处理,报错
        //子问题目标函数的一次项
        Vector *Gxk = vector_alloc(cons->dim);
        matrix_mutiply_vector(G, x_k, Gxk);
        vector_print(Gxk);
        Vector *Gxk_c = vector_alloc(cons->dim);
        vector_add_vector(Gxk, c, Gxk_c);
        Vector *p = vector_alloc(cons->dim);
        vector_print(Gxk_c);
        //子问题约束矩阵
        Matrix *sub_A = matrix_alloc(cons->dim, cons->dim);
        Vector *sub_b = vector_alloc(cons->dim);
        constrains_subconstrains(cons, W_k, sub_A, sub_b);
        vector_fill_const(sub_b, 0);
        printf("constrains matrix of subproblem when k = %d\n", k);
        matrix_print(sub_A);
        printf("constrains vector of subproblem when k = %d\n", k);
        vector_print(sub_b);
        //计算子问题得到p
        optimize_qp_linear_constraints(G, Gxk_c, sub_A, sub_b, p, y);
        printf("p in subproblem\n");
        vector_print(p);

        if (1)
        {
            //* 计算拉格朗日系数 lambda i
            // sum_{i\in W}{a_i \lambdai = g = Gx +c}

            // 子线性方程组:
            Matrix *sub_Ai = matrix_alloc(index_set_size(W_k), index_set_size(W_k));
            matrix_submatrix_by_rowindex_set(cons->A, W_k, sub_Ai);
            Vector *lambda = vector_alloc(index_set_size(W_k));
            Vector *sub_bi = vector_alloc(index_set_size(W_k));
            linear_equation_gauss_sidel(sub_Ai, sub_bi, 0.001, x0, x_k);
            if (1) //若lambda i >0 (激活不等式约束集)
            {
                ; //停止迭代
                break;
            }
            else
            {
                ; //删除某个约束
            }
        }
        else
        {
            //计算alphak
            //更新xk
            if (1) //若不满足某些约束
            {
                ; //将不满足的约束添加进工作集
            }
            else
            {

                ; //约束集不变
            }
        }
        k++;
    }

    // free workspace
    ;
    return 0;
}