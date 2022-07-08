#include "vector.h"
#include "matrix.h"
#include "util.h"
#include "math.h"
#define MAXITER 500

void __one_line(const Matrix *mat, const Vector *b, Vector *x)
{
    Vector *alpha = vector_alloc(mat->row_size);
    for (int i = 0; i < b->size; i++)
    {
        alpha->entry[i] = matrix_get_value(mat, i, 0);
    }
    double temp = vector_inner_product(alpha, b) / (vector_inner_product(alpha, alpha));
    vector_free(alpha);
    x->entry[0] = temp;
    return;
}

void linear_equation_jacobi(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x)
{
    // see: https://blog.csdn.net/Reborn_Lee/article/details/80959509
    if (!(mat->row_size == x->size AND x->size == b->size AND x->size == b->size))
    {
        terminate("ERROR :linear_equation_jacobi:input argument size does't fit.");
    }
    // x = -D^{-1}(L+U)x + D^{-1}b;
    Vector *tempx = vector_alloc(x->size);
    vector_copy(x0, tempx);
    // tempx = x;
    double s;
    long iter = 0;
    while (1)
    {
        // compute
        for (int i = 0; i < b->size; i++)
        {
            s = 0.0;
            for (int j = 0; j < b->size; j++)
            {
                if (j == i)
                {
                    continue;
                }
                s += matrix_get_value(mat, i, j) * tempx->entry[j];
            }
            x->entry[i] = (b->entry[i] - s) / matrix_get_value(mat, i, i);
        }
        iter++;
        // check iter stop
        if (vector_2metric(tempx, x) <= epsilon)
        {
            break;
        }
        // if iter bigger than MAXITER
        if (iter >= MAXITER)
        {
            terminate("ERROR :linear_equation_jacobi does't convergence");
        }
        // copy tempx
        vector_copy(x, tempx);
    }
    vector_free(tempx);
    return;
}

void linear_equation_gauss_sidel(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x)
{
    // see : https://blog.csdn.net/ITGGU/article/details/122169235
    if (!(mat->row_size == x->size AND x->size == b->size AND x->size == b->size))
    {
        terminate("ERROR :linear_equation_gauss_sidel:input argument size does't fit.");
    }
    Vector *xk = vector_alloc(x->size);
    Vector *xk_1 = vector_alloc(x->size);
    vector_copy(x0, xk);   // copy x0 to xk;
    vector_copy(x0, xk_1); // copy x0 to xk;

    double s1, s2;
    int iter = 0;
    while (1)
    {
        for (int i = 0; i < xk->size; i++) // compute
        {
            s1 = 0.0;
            s2 = 0.0;
            for (int j = 0; j < i; j++)
            {
                s1 = s1 + mat->matrix_entry[i][j] * xk_1->entry[j];
            }
            for (int j = i + 1; j < xk->size; j++)
            {
                s2 = s2 + mat->matrix_entry[i][j] * xk->entry[j];
            }
            xk_1->entry[i] = (b->entry[i] - s1 - s2) / mat->matrix_entry[i][i]; // todo
        }
        iter++;
        // check iter stop
        if (vector_2metric(xk_1, xk) <= epsilon)
        {
            break;
        }

        // if iter bigger than MAXITER
        if (iter >= MAXITER)
        {
            terminate("ERROR :linear_equation_gauss_sidel doesn't convergence");
        }
        // copy tempx
        vector_copy(xk_1, xk);
    }
    vector_copy(xk_1, x);
    vector_free(xk);
    vector_free(xk_1);
}

void linear_equation_sor(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x, const double omega)
{
    // see : https://zhuanlan.zhihu.com/p/148287671
    // omega 是超松弛系数 0< omega <2
    if (!(mat->row_size == x->size AND x->size == b->size AND x->size == b->size))
    {
        terminate("ERROR :linear_equation_sor:input argument size doesn't fit.");
    }
    Vector *xk = vector_alloc(x->size);
    Vector *xk_1 = vector_alloc(x->size);
    vector_copy(x0, xk);   // copy x0 to xk;
    vector_copy(x0, xk_1); // copy x0 to xk;

    double s1, s2;
    int iter = 0;
    while (1)
    {
        for (int i = 0; i < xk->size; i++) // compute
        {
            s1 = 0.0;
            s2 = 0.0;
            for (int j = 0; j < i; j++)
            {
                s1 = s1 + mat->matrix_entry[i][j] * xk_1->entry[j];
            }
            for (int j = i + 1; j < xk->size; j++)
            {
                s2 = s2 + mat->matrix_entry[i][j] * xk->entry[j];
            }
            xk_1->entry[i] = omega * (b->entry[i] - s1 - s2) / mat->matrix_entry[i][i] + (1. - omega) * xk->entry[i];
        }
        iter++;
        // check iter stop
        if (vector_2metric(xk_1, xk) <= epsilon)
        {
            break;
        }
        // if iter bigger than MAXITER
        if (iter >= MAXITER)
        {
            terminate("ERROR :linear_equation_sor doesn't convergence");
        }
        // copy tempx
        vector_copy(xk_1, xk);
    }
    vector_copy(xk_1, x);
    vector_free(xk);
    vector_free(xk_1);
}

int _get_max_i_number(const Matrix *A, int k)
{
    // 获取主元所在行
    int laber = k;
    double maxinum = 0;
    for (int i = k; i < A->row_size; i++)
    {
        if (maxinum < fabs(A->matrix_entry[i][k]))
        {
            maxinum = fabs(A->matrix_entry[i][k]);
            laber = i;
        }
    }
    return laber;
}

void _swap_row_ij(Matrix *A, int i, int j)
{
    // 交换矩阵第i行和第j行
    double temp = 0;
    for (int k = 0; k < A->col_size; k++)
    {
        temp = A->matrix_entry[i][k];
        A->matrix_entry[i][k] = A->matrix_entry[j][k];
        A->matrix_entry[j][k] = temp;
    }
}

void linear_equation_gaussian_elimination(const Matrix *mat, const Vector *b, Vector *x) // A[N][M]
{
    // see: https://blog.csdn.net/weixin_42465397/article/details/103264328
    if (!(mat->col_size == x->size AND mat->row_size == b->size))
    {
        terminate("ERROR :linear_equation_gaussian_elimination:input argument size doesn't fit.");
    }
    // 处理一个特殊情况,矩阵行为1的情况:
    if (mat->col_size == 1)
    {
        __one_line(mat, b, x);
        return;
    }

    int laber;
    double temp;
    double sum;
    Matrix *tempA = matrix_alloc(mat->row_size, mat->col_size + 1);
    // 构建增广矩阵
    // N:row_size
    // M:col_size
    vector_fill_const(x, 0.0);
    for (int i = 0; i < tempA->row_size; i++)
    {
        for (int j = 0; j < tempA->col_size; j++)
        {
            if (j < tempA->col_size - 1)
                // tempA->matrix_entry[i][j] = mat->matrix_entry[i][j];
                matrix_set_value(tempA, i, j, matrix_get_value(mat, i, j));
            else
                // tempA->matrix_entry[i][j] = b->entry[i];
                matrix_set_value(tempA, i, j, b->entry[i]);
        }
    }
    for (int k = 0; k < mat->row_size; k++)
    {
        laber = _get_max_i_number(tempA, k);
        if (laber != k)
        {
            _swap_row_ij(tempA, laber, k);
        }

        for (int i = k + 1; i < mat->row_size; i++)
        {
            // if (tempA->matrix_entry[k][k] == 0.0)
            if (matrix_get_value(tempA, k, k) == 0.0)
                break;
            // temp = tempA->matrix_entry[i][k] / tempA->matrix_entry[k][k];
            temp = matrix_get_value(tempA, i, k) / matrix_get_value(tempA, k, k);
            for (int j = k; j < tempA->col_size; j++)
                // tempA->matrix_entry[i][j] = tempA->matrix_entry[k][j] * temp - tempA->matrix_entry[i][j];
                matrix_set_value(tempA, i, j,
                                 matrix_get_value(tempA, k, j) * temp - matrix_get_value(tempA, i, j));
        }
    }
    for (int i = tempA->row_size - 1; i >= 0; i--)
    {
        sum = 0.;
        for (int j = i + 1; j < tempA->row_size; j++)
        {
            sum = sum + tempA->matrix_entry[i][j] * x->entry[j];
        }
        x->entry[i] = (tempA->matrix_entry[i][tempA->col_size - 1] - sum) / tempA->matrix_entry[i][i];
    }
    matrix_free(tempA);
}