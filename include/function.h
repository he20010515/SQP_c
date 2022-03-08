#include "myvector.h"
#pragma once
struct NdsclaFunction
{
    unsigned int inputSize;
    double (*function)(Vector *);
};
typedef struct NdsclaFunction NdsclaFunction;
float NdsclaFunctionCall(NdsclaFunction *function, Vector *x);
NdsclaFunction *NdsclaFunctionAlloc(double (*function)(Vector *), int inputsize);
void centralGrad(NdsclaFunction *function, double h, Vector *x0, Vector *grad);
