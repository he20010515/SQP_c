#include "vector.h"
#include <stdlib.h>
#include "function.h"
#include "matrix.h"
double NdsclaFunctionCall(NdsclaFunction *function, Vector *x)
{
    return function->function(x);
}

NdsclaFunction *NdsclaFunctionAlloc(double (*function)(Vector *), int inputsize)
{
    NdsclaFunction *f = (NdsclaFunction *)malloc(sizeof(NdsclaFunction));
    f->function = function;
    f->inputSize = inputsize;
}

void centralGrad(NdsclaFunction *function, double h, Vector *x0, Vector *grad)
{

    Vector *temp = VectorAlloc(x0->size);
    VectorCopy(x0, temp);
    double f_add_h, f_sub_h;
    for (size_t i = 0; i < function->inputSize; i++)
    {
        temp->entry[i] += h;
        f_add_h = NdsclaFunctionCall(function, temp);
        temp->entry[i] -= 2 * h;
        f_sub_h = NdsclaFunctionCall(function, temp);
        temp->entry[i] += h;
        grad->entry[i] = (f_add_h - f_sub_h) / (2. * h);
    }
    VectorFree(temp);
}

void centralHession(NdsclaFunction *function, double h, Vector *x0, Matrix *hession)
{
    // 计算Hession矩阵 H[i][j] = \frac{\partial f}{\partial x_i \partial x_j}
    if (hession->col_size != hession->row_size)
        terminate("ERROR hession matrix must have the same row and col size");
    else if (hession->col_size != function->inputSize)
        terminate("ERROT function inputsize must same as hession size");
    for (size_t i = 0; i < hession->row_size; i++)
    {
        for (size_t j = 0; j < hession->row_size; j++)
        {

            hession->matrix_entry[i][j] = 0;
        }
    }
}