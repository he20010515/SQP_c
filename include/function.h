#include "vector.h"
#pragma once
struct NdsclaFunction
{
    unsigned int inputSize;
    double (*function)(Vector *);
};
typedef struct NdsclaFunction NdsclaFunction;
double NdsclaFunctionCall(NdsclaFunction *function, Vector *x);
NdsclaFunction *NdsclaFunctionAlloc(double (*function)(Vector *), int inputsize);
void centralGrad(NdsclaFunction *function, double h, Vector *x0, Vector *grad);
