#include "matrix.h"
#include "vector.h"
#include "function.h"
#include "util.h"
#include "linear_equations.h"
int optimize_qp(const Matrix *H, const Vector *c, const Matrix *A, const Vector *b, Vector *x_star, double *maxy)
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
        terminate("ERROR: optimize_qp:input argument does NOT fit, please check input");
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
    //* compute linear function #TODO 这里用一般的迭代法不收敛,目前先用LU分解验证
    // linear_equation_sor(HAAT0, cb, 0.0001, xlambda0, xlambda, 0.5);
    Matrix *inv = matrix_alloc(2 * n, 2 * n);
    matrix_inverse(HAAT0, inv);
    matrix_mutiply_vector(HAAT0, cb, xlambda);
    matrix_free(inv);

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