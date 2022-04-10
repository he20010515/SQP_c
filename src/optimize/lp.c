#include "vector.h"
#include "matrix.h"
#include "lp.h"
#include "util.h"
#include "math.h"
#include "simplex.h"
#include "elog.h"

#define LOG_TAG "lp"
/*
求解线性规划(LP)问题的程序,主要用于计算QP问题中的初始值
*/

//线性约束
LinearConstraints *linearconstraints_alloc(int dim, int size, int e, int i, Matrix *A, Vector *b)
{
    LinearConstraints *con = (LinearConstraints *)malloc(sizeof(LinearConstraints));
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
    return con;
}

void *linearconstraints_verification(const LinearConstraints *con, const Vector *x, Index_set *set)
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

void linearconstraints_free(LinearConstraints *con, int recursion)
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

void linearconstrains_subconstrains(const LinearConstraints *con, const Index_set *set, Matrix *A, Vector *b)
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

int optimize_lp(const LinearConstraints *con, const Vector *c, Vector *x0, int maxiter, double tol, int bland)
{
    // Trans Problem From:
    // min c^T x
    // s.t. Ax == b (i < e    )
    //      Ax >= b (e <=i < m)
    // TO:
    // min c^T (x-x') + 0^T *xs
    // s.t.   A(x-x')      == b
    //        A(x-x') - xs == b
    // Use simplex method to solve this problem

    // matrixA' :
    // [        [ 0 ]]
    // [A  -A   [ E ]]
    // [        [ E ]]
    // X" = [x0',x1',...,xn-1',x0",x1",...,xn-1",xs1,xs2,..,xsm];
    int n = con->dim;  //维度
    int m = con->size; //原约束大小
    int k = con->e;    //等式数量
    Matrix *mat = matrix_alloc(m, 2 * n + m - k);
    for (int i = 0; i < mat->row_size; i++)
    {
        for (int j = 0; j < mat->col_size; j++)
        {
            if (j < 2 * n)
            {
                // A -A
                if (j < n)
                    mat->matrix_entry[i][j] = con->A->matrix_entry[i][j];
                if (n <= j AND j < 2 * n)
                    mat->matrix_entry[i][j] = -con->A->matrix_entry[i][j - n];
            }
            else // j>=2*n
            {
                if (i < k)
                    mat->matrix_entry[i][j] = 0.0;
                else
                {
                    if (i - k == j - 2 * n)
                        mat->matrix_entry[i][j] = -1;
                    else
                        mat->matrix_entry[i][j] = 0;
                }
            }
        }
    }
    Vector *b = vector_alloc(con->b->size);
    for (int i = 0; i < mat->row_size; i++)
    {
        if (con->b->entry[i] < 0)
        {
            b->entry[i] = -con->b->entry[i];
            for (int j = 0; j < mat->col_size; j++)
                mat->matrix_entry[i][j] *= -1;
        }
        else
            b->entry[i] = con->b->entry[i];
    }
    Vector *tempc = vector_alloc(mat->col_size);
    for (int i = 0; i < mat->col_size; i++)
    {
        if (i < 2 * n)
        {
            if (i < n)
                tempc->entry[i] = c->entry[i];
            else
                tempc->entry[i] = -c->entry[i - con->dim];
        }
        else
            tempc->entry[i] = 0.;
    }
    Vector *temp_x0 = vector_alloc(2 * n + m - k);
    // matrix_print(mat);
    // vector_print(tempc);
    // vector_print(b);
    int flag = _linprog_simplex(tempc, mat, b, maxiter, tol, bland, temp_x0);
    if (flag == 0)
    {
        // case optional solution found
        for (int i = 0; i < con->dim; i++)
            x0->entry[i] = temp_x0->entry[i] - temp_x0->entry[i + con->dim];
    }
    matrix_free(mat);
    vector_free(b);
    vector_free(tempc);
    vector_free(temp_x0);
    return flag;
}