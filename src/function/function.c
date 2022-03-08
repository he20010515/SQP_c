#include "vector.h"
#include <stdlib.h>
#include "function.h"
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
