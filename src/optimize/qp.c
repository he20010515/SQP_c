#include "matrix.h"
#include "vector.h"
#include "function.h"
#include "util.h"
#include "linear_equations.h"
#include "index_set.h"
#include "qp.h"
#include "math.h"

#define min(a, b) (((a) < (b)) ? (a) : (b))

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
    if (!(con->size == set->index_range AND A->col_size == con->dim AND b->size == index_set_size(set) AND A->row_size == index_set_size(set)))
    { // 约束的大小等于索引集合的范围;约束矩阵的列数等于问题维度;约束列向量大小等于约束数量;矩阵行数等于工作集大小
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
    Index_set *index_set_I = index_set_alloc(m);
    for (int i = cons->e; i < m; i++)
    {
        index_set_append(index_set_I, i);
    }
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
        printf("========iter k = %d=========\n", k);
        Vector *p = vector_alloc(cons->dim);
        __qp_compute_subproblem(W_k, cons, G, c, x_k, p, y);
        if (double_equal(vector_2norm(p), 0.0)) // if p_k = 0
        {
            //*计算lambda
            Vector *lambda = vector_alloc(cons->size);
            Vector *subsublambda = __qp_compute_lambda(W_k, cons, G, c, x_k, index_set_I, lambda);
            if (vector_any_bigger_equal_than_const(subsublambda, 0)) //若lambda i >0 (激活不等式约束集) (\any i \in Wk)
            {
                vector_free(subsublambda);
                printf("case: iter done\n");
                vector_copy(x_k, x_star);
                return 0;
            }
            else
            {
                printf("case: remove con\n");
                int j = vector_argmin(lambda);
                index_set_remove(W_k, j);
                vector_copy(x_k, x_k_1);
            }
        }
        else
        {
            //*计算alphak
            //更新xk
            Vector *alphas = vector_alloc(cons->size);
            double alphak = __qp_compute_alphak(W_k, cons, x_k, p, alphas);
            Vector *alphapk = vector_multiply_const(p, alphak, 1); // copy = 1
            int j = vector_argmin(alphas);
            vector_add_vector(x_k, alphapk, x_k_1);
            vector_free(alphapk);
            vector_free(alphas);
            if (alphak < 1.) //若不满足某些约束
            {
                index_set_append(W_k, j); //将不满足的约束添加进工作集
                printf("Index set append%d\n", j);
            }
            else
            {
                ; //约束集不变
            }
        }
        printf("iter%d done;xk = \n", k);
        vector_print(x_k);

        k++;
        vector_copy(x_k_1, x_k);
    }

    // free workspace
    ;
    return 0;
}

void __qp_compute_subproblem(const Index_set *W_k, const Constraints *cons, const Matrix *G, const Vector *c, Vector *xk, Vector *p, double *y)
{
    int sizeofw = index_set_size(W_k); // TODO 将等式约束加入指标集  DONE
    for (int i = 0; i < cons->e; i++)  //确保添加进去等式约束
        index_set_append(W_k, i);
    if (sizeofw > cons->dim)
        terminate("ERROR size of wk too big"); //! 如果指标集的大小比输入维度还大的话,这种问题我们目前还没法处理,报错
    //子问题目标函数的一次项
    Vector *Gxk = vector_alloc(cons->dim);
    matrix_mutiply_vector(G, xk, Gxk);
    Vector *Gxk_c = vector_alloc(cons->dim);
    vector_add_vector(Gxk, c, Gxk_c);
    //子问题约束矩阵
    Matrix *sub_A = matrix_alloc(index_set_size(W_k), cons->dim); // 约束个行,dim个列
    Vector *sub_b = vector_alloc(index_set_size(W_k));
    constrains_subconstrains(cons, W_k, sub_A, sub_b);
    vector_fill_const(sub_b, 0);
    //计算子问题得到p
    optimize_qp_linear_constraints(G, Gxk_c, sub_A, sub_b, p, y);
    return;
}
Vector *__qp_compute_lambda(const Index_set *W_k, const Constraints *cons, const Matrix *G, const Vector *c, Vector *xk, Index_set *index_set_I, Vector *lambda)
{
    //* 计算拉格朗日系数 lambda i
    // sum_{i\in W}{a_i \lambdai = g = Gx +c}
    // 子线性方程组:
    int m = cons->dim;
    Vector *Gxk = vector_alloc(cons->dim);
    matrix_mutiply_vector(G, xk, Gxk);
    Vector *Gxk_c = vector_alloc(cons->dim);
    vector_add_vector(Gxk, c, Gxk_c);
    Matrix *sub_Ai = matrix_alloc(index_set_size(W_k), cons->dim);
    Matrix *sub_AiT = matrix_alloc(cons->dim, index_set_size(W_k));
    matrix_submatrix_by_rowindex_set(cons->A, W_k, sub_Ai);
    matrix_transpose(sub_Ai, sub_AiT);
    Vector *sub_lambda = vector_alloc(index_set_size(W_k));           // wk大小的lambda
    linear_equation_gaussian_elimination(sub_AiT, Gxk_c, sub_lambda); //小的lambda
    //将lambda 放大到根activeset一样大
    int temp = 0;
    for (int i = 0; i < W_k->index_range; i++)
    {
        if (index_set_is_in(W_k, i))
        {
            lambda->entry[i] = sub_lambda->entry[temp];
            temp++;
        }
        else
            lambda->entry[i] = 0.0;
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
    return subsublambda;
}
double __qp_compute_alphak(const Index_set *W_k, const Constraints *cons, const *x_k, const Vector *p, Vector *alphas)
{
    //计算alphak
    double alphak = 1.0;
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