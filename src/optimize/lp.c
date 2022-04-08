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

int __judge(const Matrix *mat);
int __trans(Matrix *mat, int *vect);

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

double optimize_lp_standard_type(const Vector *c, const Vector *b, const Matrix *A, const int *init_base, Vector *xstar)
{
    // 求解线性规划标准形式,需要初始基集合,初始基集合必须为单位阵
    // see: https://blog.csdn.net/qq_47723068/article/details/109537450
    // Solve Problem:
    //  max z = c^Tx
    //  Ax = b
    //  x_i >=0  i = 1,2,...,n
    // check input
    if (!(c->size == A->col_size AND A->row_size == b->size))
        terminate("optimize_lp_standard_type: input don't fit");

    // 申请单纯形表矩阵
    Matrix *mat = matrix_alloc(A->row_size + 1, A->col_size + 1);
    // init
    for (int i = 0; i < A->row_size; i++)
        for (int j = 0; j < A->col_size; j++)
            mat->matrix_entry[i][j] = A->matrix_entry[i][j];
    for (int i = 0; i < b->size; i++)
        mat->matrix_entry[i][mat->col_size - 1] = b->entry[i];
    for (int i = 0; i < c->size; i++)
        mat->matrix_entry[mat->row_size - 1][i] = c->entry[i];
    mat->matrix_entry[mat->row_size - 1][mat->col_size - 1] = 0;

    int *vect = (int *)malloc(sizeof(int) * A->row_size);
    // TODO initVect
    for (int i = 0; i < A->row_size; i++)
        vect[i] = init_base[i];
    matrix_print(mat);
    while (__judge(mat))
    {
        __trans(mat, vect);
        matrix_print(mat);
    }
    // return answer:
    for (int i = 0; i < A->col_size; i++)
    {
        int flag = -1;
        for (int j = 0; j < A->row_size; j++)
        {
            if (vect[j] == i)
            {
                flag = j;
                break;
            }
            else

                flag = -1;
        }
        // if i in vect:
        if (flag != -1)
            xstar->entry[i] = mat->matrix_entry[flag][mat->col_size - 1];
        // else
        else
            xstar->entry[i] = 0.0;
    }
    matrix_free(mat);
    return vector_inner_product(c, xstar);
}

int __judge(const Matrix *mat)
{
    double m = mat->matrix_entry[mat->row_size - 1][0];
    for (int i = 0; i < mat->col_size - 1; i++)
        m = MAX(m, mat->matrix_entry[mat->row_size - 1][i]);
    if (m <= 0.0)
        return FALSE;
    else
        return TRUE;
}

int __trans(Matrix *mat, int *vect)
{
    Vector *temp = vector_alloc(mat->col_size - 1);
    for (int i = 0; i < mat->col_size - 1; i++)
        temp->entry[i] = mat->matrix_entry[mat->row_size - 1][i];
    int in_base = vector_argmax(temp); // 入基变量

    vector_free(temp);
    Vector *temp_2 = vector_alloc(mat->row_size - 1);
    vector_fill_const(temp_2, NAN);
    for (int i = 0; i < mat->row_size - 1; i++) // 不遍历最后一行
        if (mat->matrix_entry[i][in_base] > 0.0)
            temp_2->entry[i] = mat->matrix_entry[i][mat->col_size - 1] / mat->matrix_entry[i][in_base];
    int out_base = vector_argmin(temp_2); //出基变量的角标;
    vector_free(temp_2);

    double k = (mat->matrix_entry[out_base][in_base]);
    for (int i = 0; i < mat->col_size; i++) // out_base[pivot][index]
        mat->matrix_entry[out_base][i] = ((mat->matrix_entry[out_base][i]) / k);
    for (int i = 0; i < mat->row_size; i++)
    {
        if (i != out_base)
        {
            // matrix的第i行等于 matrix的第i行减去 matrix[i][index] *matrix[prvot]
            k = mat->matrix_entry[i][in_base];
            for (int j = 0; j < mat->col_size; j++)
                mat->matrix_entry[i][j] = mat->matrix_entry[i][j] - k * mat->matrix_entry[out_base][j];
        }
    }
    vect[out_base] = in_base;
}

int optimize_lp(const LinearConstraints *con, const Vector *c, Vector *x0)
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
    Simplex *problem = simplex_alloc(tempc, NULL, NULL, mat, b);
    Vector *temp_x0 = vector_alloc(2 * n + m - k);

    // matrix_print(mat);
    // vector_print(tempc);
    // vector_print(b);

    int flag = simplex_main(problem, temp_x0);

    if (flag == 1)
    {
        // case optional solution found
        for (int i = 0; i < con->dim; i++)
            x0->entry[i] = temp_x0->entry[i] - temp_x0->entry[i + con->dim];
    }
    else
    {
        // case infity solution
        for (int i = 0; i < con->dim; i++)
            x0->entry[i] = INFINITY;
    }

    simplex_free(problem);
    matrix_free(mat);
    vector_free(b);
    vector_free(tempc);
    vector_free(temp_x0);
    return flag;
}