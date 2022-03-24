#include "vector.h"
#include <stdlib.h>
#include "function.h"
#include "matrix.h"
#include "elog.h"
double ndscla_function_call(NdsclaFunction *function, Vector *x)
{
    return function->function(x);
}

NdsclaFunction *ndscla_function_alloc(double (*function)(Vector *), int inputsize)
{
    NdsclaFunction *f = (NdsclaFunction *)malloc(sizeof(NdsclaFunction));
    f->function = function;
    f->inputSize = inputsize;
    return f;
}

void ndscla_central_grad(NdsclaFunction *function, double h, Vector *x0, Vector *grad)
{

    Vector *temp = vector_alloc(x0->size);
    vector_copy(x0, temp);
    double f_add_h, f_sub_h;
    for (size_t i = 0; i < function->inputSize; i++)
    {
        temp->entry[i] += h;
        f_add_h = ndscla_function_call(function, temp);
        temp->entry[i] -= 2 * h;
        f_sub_h = ndscla_function_call(function, temp);
        temp->entry[i] += h;
        grad->entry[i] = (f_add_h - f_sub_h) / (2. * h);
    }
    vector_free(temp);
}

void ndscla_central_hession(NdsclaFunction *function, double h, Vector *x0, Matrix *hession)
{
    // 计算Hession矩阵 H[i][j] = \frac{\partial f}{\partial x_i \partial x_j}
    if (hession->col_size != hession->row_size)
        terminate("ERROR hession matrix must have the same row and col size");
    else if (hession->col_size != function->inputSize)
        terminate("ERROT function inputsize must same as hession size");

    // begincompute
    double f_ai_aj, f_mi_aj, f_ai_mj, f_mi_mj;
    Vector *tempv = vector_alloc(function->inputSize);
    vector_copy(x0, tempv);
    for (size_t i = 0; i < hession->row_size; i++)
    {
        for (size_t j = 0; j < hession->row_size; j++)
        {
            tempv->entry[i] += h; // x+hi
            tempv->entry[j] += h; // x+hi+hj
            f_ai_aj = ndscla_function_call(function, tempv);
            tempv->entry[i] -= 2 * h; // x-hi+hj
            f_mi_aj = ndscla_function_call(function, tempv);
            tempv->entry[j] -= 2 * h; // x-hi-hj
            f_mi_mj = ndscla_function_call(function, tempv);
            tempv->entry[i] += 2 * h; // x+hi-hj
            f_ai_mj = ndscla_function_call(function, tempv);
            tempv->entry[i] -= h;
            tempv->entry[j] += h;

            hession->matrix_entry[i][j] = (f_ai_aj - f_ai_mj - f_mi_aj + f_mi_mj) / (4 * h * h);
        }
    }
    vector_free(tempv);
    return;
}

void ndsclaFunction_free(NdsclaFunction *function)
{
    free(function);
}

NdVectorfunction *ndVectorfunction_alloc(void (*function)(Vector *, Vector *), int intputdim, int outputdim)
{
    NdVectorfunction *f = (NdVectorfunction *)malloc(sizeof(NdVectorfunction));
    f->function = function;
    f->inputdim = intputdim;
    f->outputdim = outputdim;
}

void ndVectorfunction_free(NdVectorfunction *function)
{
    free(function);
}

void ndVectorfunction_call(NdVectorfunction *function, Vector *input, Vector *output)
{
    if (!(input->size == function->inputdim AND output->size == function->outputdim))
    {
        terminate("shape not fit");
    }
    function->function(input, output);
    return;
}

void ndVectorfunction_jacobian(const NdVectorfunction *function, const Vector *x0, double h, Matrix *jacobian)
{
    Vector *y = vector_alloc(function->outputdim);
    Vector *x = vector_alloc(function->inputdim);
    double yi_add_h, yi_sub_h;
    vector_copy(x0, x);
    if (!(jacobian->col_size == function->inputdim AND jacobian->row_size == function->outputdim AND x0->size == function->inputdim))
    {
        terminate("ndVectorfunction_jacobian size don't fit");
    }
    for (int i = 0; i < jacobian->row_size; i++)
    {
        for (int j = 0; j < jacobian->col_size; j++)
        {
            // \partical yi/ \partical xj;
            ;
            x->entry[j] += h;
            ndVectorfunction_call(function, x, y);
            yi_add_h = y->entry[i];
            x->entry[j] -= 2 * h;
            ndVectorfunction_call(function, x, y);
            yi_sub_h = y->entry[i];
            x->entry[j] += h;

            jacobian->matrix_entry[i][j] = (yi_add_h - yi_sub_h) / (2. * h);
        }
    }
    vector_free(y);
    vector_free(x);
}
