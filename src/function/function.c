#include "vector.h"
#include <stdlib.h>
#include "function.h"
#include "matrix.h"
#include "elog.h"
#include <math.h>
#include <omp.h>
#include "sparse_matrix.h"

double ndscla_function_call(const NdsclaFunction *function, Vector *x)
{
#ifndef FUNCTION_RECORED
    NdsclaFunction *f = (NdsclaFunction *)function;
    f->real_call_num++;
    return function->function(x);
#else
    NdsclaFunction *f = (NdsclaFunction *)function;
    double *y = (double *)HashTable_get(f->table, x);
    if (y != NULL)
    {
        f->record_call_num++;
        return *y;
    }
    else
    {
        f->real_call_num++;
        y = (double *)malloc(sizeof(double));
        *y = function->function(x);
        Vector *copy_x = vector_alloc(x->size);
        vector_copy(x, copy_x);
        Pointer_buffer_insert(function->x_buffer, copy_x);
        Pointer_buffer_insert(function->y_buffer, y);
        HashTable_insert(function->table, copy_x, y);
        return *y;
    }

#endif
}

#ifdef FUNCTION_RECORED
unsigned int vector_hash(void *vector)
{
    Vector *v = (Vector *)vector;
    unsigned int hash = 114514;
    char *key = (char *)v->entry;
    for (int i = 0; i < v->size * sizeof(double); key++, i++)
    {
        hash ^= ((hash << 5) + (*key) + (hash >> 2));
    }
    return hash;
}
int vector_compare(void *_v1, void *_v2)
{
    Vector *v1 = (Vector *)_v1, *v2 = (Vector *)_v2;
    double t = 0.0;
    for (int i = 0; i < v1->size; i++)
    {
        t += (v1->entry[i] - v2->entry[i]) * (v1->entry[i] - v2->entry[i]);
    }
    t = sqrt(t);
    return t <= 1e-12;
}

#endif

NdsclaFunction *ndscla_function_alloc(double (*function)(Vector *), int inputsize)
{
    NdsclaFunction *f = (NdsclaFunction *)malloc(sizeof(NdsclaFunction));
    f->real_call_num = 0;
    f->record_call_num = 0;
    f->function = function;
    f->inputSize = inputsize;
    f->table = NULL;
    f->x_buffer = NULL;
    f->y_buffer = NULL;
#ifdef FUNCTION_RECORED
    f->table = HashTable_alloc(vector_hash, vector_compare);
    f->x_buffer = Pointer_buffer_alloc();
    f->y_buffer = Pointer_buffer_alloc();
#endif
    return f;
}

void ndscla_function_free(NdsclaFunction *f)
{
    free(f);
#ifdef FUNCTION_RECORED
    HashTable_free(f->table);
    Pointer_buffer_free(f->x_buffer);
    Pointer_buffer_free(f->y_buffer);
#endif
}

void ndscla_forward_grad(const NdsclaFunction *function, double h, const Vector *x0, Vector *grad)
{

    Vector *temp = vector_alloc(x0->size);
    vector_copy(x0, temp);
    int i = 0;
    double f_add_h, f;
    //#pragma omp parallel for num_threads(10) default(none) shared(function, x0, grad, temp) private(i, h, f_add_h, f)
    for (i = 0; i < function->inputSize; i++)
    {
        h = 2 * (1 + fabs(temp->entry[i])) * sqrtl(exp2l(log2l(fabs(temp->entry[i]) + 2e-10) - 55.0));
        temp->entry[i] += h;
        f_add_h = ndscla_function_call(function, temp);
        temp->entry[i] -= h;
        f = ndscla_function_call(function, temp);
        grad->entry[i] = (f_add_h - f) / (h);
        // printf("compute from thread %3d \n", omp_get_thread_num());
    }
    vector_free(temp);
}

void ndscla_central_hession(const NdsclaFunction *function, double h, const Vector *x0, Matrix *hession)
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
    for (int i = 0; i < hession->row_size; i++)
    {
        for (int j = 0; j < hession->row_size; j++)
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

NdVectorfunction *ndVectorfunction_alloc(void (*function)(const Vector *, Vector *), int intputdim, int outputdim)
{
    NdVectorfunction *f = (NdVectorfunction *)malloc(sizeof(NdVectorfunction));
    f->function = function;
    f->inputdim = intputdim;
    f->outputdim = outputdim;
    return f;
}

void ndVectorfunction_free(NdVectorfunction *function)
{
    free(function);
}

void ndVectorfunction_call(const NdVectorfunction *function, const Vector *input, Vector *output)
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