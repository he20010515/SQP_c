#include "vector.h"
#include "matrix.h"
#include "lp.h"
#include "util.h"
#include "math.h"
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

//求解线性规划标准形式, 不需要初始解
void optimize_lp_2stage(const Vector *c, const Vector *b, const Matrix *A, Vector *xstar)
{
    // see: https://baike.baidu.com/item/%E4%B8%A4%E9%98%B6%E6%AE%B5%E6%B3%95/9302919
    // Solve Problem:
    //  max z = c^Tx
    //  Ax = b
    //  x_i >=0  i = 1,2,...,n
    int m = A->row_size;
    int n = A->col_size;
    // step1 : 添加人工变量
    //  min w = 0^Tx + 1^T xr
    //  Ax + Exr = b
    //  x_i >=0  i = 1,2,...,n
    Vector *xxr = vector_alloc(m + n);   //原维度+约束数量
    Matrix *AA = matrix_alloc(m, m + n); // step1约束矩阵
    const Vector *bb = b;
    Vector *cc = vector_alloc(m + n);
    for (int i = 0; i < m + n; i++) // fill cc
    {
        if (i < n)
            cc->entry[i] = 0.0;
        else
            cc->entry[i] = 1.0;
    }
    for (int i = 0; i < m; i++) // fill AA
    {
        for (int j = 0; j < m + n; j++)
        {
            if (j < m)
            {
                AA->matrix_entry[i][j] = A->matrix_entry[i][j];
            }
            else
            {
                if (i == j - n)
                    AA->matrix_entry[i][j] = 1.0;
                else
                    AA->matrix_entry[i][j] = 0.0;
            }
        }
    }
    int *vect = (int *)malloc(sizeof(int) * m); // 初始基
    for (int i = 0; i < m; i++)
        vect[i] = n + i;
    vector_print(cc);
    matrix_print(AA);
    optimize_simplex_method(cc, bb, AA, xxr, vect);
    vector_print(xxr);
    // step2 计算解
    free(vect);
}

int __find_swap_in_base(int n, Vector *sigma);
int __find_swap_out_base(int a, Matrix *A, Vector *seta, Vector *b, int *num, Vector *CB, Vector *C);
int __iterate(int p, int q, Matrix *A, Vector *b, Vector *sigma);

void optimize_simplex_method(const Vector *_c, const Vector *_b, const Matrix *_A, Vector *xstar, int *init_base)
{
    // see :https://blog.csdn.net/weixin_47000625/article/details/109341360
    // Solve Problem
    //  min   c^Tx
    //  s.t.  Ax = b
    if (!(_c->size == _A->col_size AND _A->row_size == _b->size))
        terminate("optimize_lp_standard_type: input don't fit");
    int m = _A->row_size;
    int n = _A->col_size;
    // alloc workspace
    Matrix *A = matrix_alloc(m, n);  //约束矩阵
    Vector *C = vector_alloc(n);     //目标函数价值系数
    Vector *b = vector_alloc(m);     //资源系数
    Vector *CB = vector_alloc(m);    //基变量系数矩阵
    Vector *seta = vector_alloc(m);  //存放出基和入基的变化情况
    Vector *sigma = vector_alloc(n); //检验数
    Vector *x = vector_alloc(n);     //决策变量
    double z = 0;                    //目标函数值
    int *num = (int *)malloc(m * sizeof(int));
    // init variable
    matrix_copy((Matrix *)_A, A);
    vector_copy(_c, C);
    vector_copy(_b, b);
    vector_fill_const(x, 0.0);
    for (int i = 0; i < m; i++)
        num[i] = init_base[i];
    vector_copy(C, sigma);
    for (int i = 0; i < m; i++)
        CB->entry[i] = C->entry[num[i]];
    int p, q; // q是换入变量
    while (TRUE)
    {
        q = __find_swap_in_base(n, sigma);
        if (q == -1)
        {
            for (int j = 0; j < m; j++)
                x->entry[num[j]] = b->entry[j];
            vector_copy(x, xstar);
            z = vector_inner_product(x, C);
            printf("z = === %lf", z);
            break;
        }
        p = __find_swap_out_base(q, A, seta, b, num, CB, C);
        __iterate(p, q, A, b, sigma);
    }
}

int __find_swap_in_base(int n, Vector *sigma)
{
    int i, k = 0;
    int flag = 1;
    double min = 0.0;
    for (i = 0; i < n; i++)
    {
        if (sigma->entry[i] < 0)
        {
            flag = 0;
            break;
        }
    }
    if (flag == 1)
        return -1;
    for (i = 0; i < n; i++)
    {
        if (sigma->entry[i] < min)
        {
            min = sigma->entry[i];
            k = i;
        }
    }
    return k;
}

int __find_swap_out_base(int a, Matrix *A, Vector *seta, Vector *b, int *num, Vector *CB, Vector *C)
{
    const double M = 1e10;
    int i, j;
    int m = A->row_size;
    int flag = 1;
    int k = a;              // a入基变量
    for (i = 0; i < m; i++) //如果某个数小于0的检验数，对应的列向量中所有元素≤0 该问题有无界解
    {
        if (A->matrix_entry[i][k] > 0)
        {
            flag = 0;
            break;
        }
    }
    if (flag == 1)
    {
        printf("该线性规划问题有无界解、无最优解！\n");
        return -1;
    }
    for (i = 0; i < m; i++)
    {
        if (A->matrix_entry[i][k] > 0)
            seta->entry[i] = b->entry[i] / A->matrix_entry[i][k];
        else
            seta->entry[i] = M; //当A[i][k]≤0的时候是不需要考虑的
                                //然而是根据比值最小原则整的 所以给对应的卡一个很大的M
    }
    double min = M;
    for (i = 0; i < m; i++) //得到换出变量
    {
        if (min >= seta->entry[i])
        {
            min = seta->entry[i];
            j = i;
        }
    }
    num[j] = k;
    CB->entry[j] = C->entry[k];
    return j;
}

int __iterate(int p, int q, Matrix *A, Vector *b, Vector *sigma)
{
    int i, j, r, c; // row行 column列（r,l)就是转轴元
    r = p;          //行号		p是出基变量
    c = q;          //列号		q是入基变量
    int n = A->col_size;
    int m = A->row_size;
    double temp1 = A->matrix_entry[r][c];
    double temp2, temp3;
    //标准化该行
    b->entry[r] /= temp1;
    for (j = 0; j < n; j++)
        A->matrix_entry[r][j] /= temp1;
    for (i = 0; i < m; i++) //标准化其他行
    {
        if (i != r)
            if (!double_equal(A->matrix_entry[i][c], 0.0))
            {
                temp2 = A->matrix_entry[i][c];
                b->entry[i] -= temp2 * b->entry[r]; // b[r]转轴元对应b
                for (j = 0; j < n; j++)
                    A->matrix_entry[i][j] -= temp2 * A->matrix_entry[r][j]; // A[r][j]是转轴元对应行
            }
    }
    //σ的迭代计算
    temp3 = sigma->entry[c];
    for (i = 0; i < n; i++)
    {
        sigma->entry[i] -= A->matrix_entry[r][i] * temp3;
    }
}