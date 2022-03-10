#include "vector.h"
#include "matrix.h"
#include "util.h"
#define MAXITER 500
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
                s += mat->matrix_entry[i][j] * tempx->entry[j];
            }
            x->entry[i] = (b->entry[i] - s) / mat->matrix_entry[i][i];
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