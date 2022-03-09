#include "vector.h"
#include "matrix.h"
#include "util.h"
#define MAXITER 10
void linear_equation_jacobi(const Matrix *mat, const Vector *b, const double epsilon, const Vector *x0, Vector *x)
{
    // see: https://blog.csdn.net/Reborn_Lee/article/details/80959509
    if (!(mat->row_size == x->size AND x->size == b->size AND x->size == b->size))
    {
        terminate("ERROR :linear_equation_jacobi:input argument size don't fit.");
    }
    // x = -D^{-1}(L+U)x + D^{-1}b;
    Vector *tempx = vector_alloc(x->size);
    vector_copy(x0, tempx);
    // tempx = x;
    double s;
    long iter = 0;
    while (1)
    {
        vector_print(x);
        // compute
        for (int i = 0; i < b->size; i++)
        {
            s = 0.0;
            for (int j = 0; j < b->size; j++)
            {
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
            terminate("ERROR :linear_equation_jacobi don't convergence");
        }
        // copy tempx
        vector_copy(x, tempx);
    }
    vector_free(tempx);
    return;
}