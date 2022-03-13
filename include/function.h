#include "vector.h"
#include "matrix.h"
#pragma once
struct NdsclaFunction
{
    unsigned int inputSize;
    double (*function)(Vector *);
};
typedef struct NdsclaFunction NdsclaFunction;
double ndscla_function_call(NdsclaFunction *function, Vector *x);
NdsclaFunction *ndscla_function_alloc(double (*function)(Vector *), int inputsize);
void ndscla_central_grad(NdsclaFunction *function, double h, Vector *x0, Vector *grad);
void ndscla_central_hession(NdsclaFunction *function, double h, Vector *x0, Matrix *hession);