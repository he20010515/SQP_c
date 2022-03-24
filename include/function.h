#include "vector.h"
#include "matrix.h"
#pragma once
struct NdsclaFunction
{
    unsigned int inputSize;
    double (*function)(Vector *);
};
typedef struct NdsclaFunction NdsclaFunction;
NdsclaFunction *ndscla_function_alloc(double (*function)(Vector *), int inputsize);
double ndscla_function_call(NdsclaFunction *function, Vector *x);
void ndscla_central_grad(NdsclaFunction *function, double h, Vector *x0, Vector *grad);
void ndscla_central_hession(NdsclaFunction *function, double h, Vector *x0, Matrix *hession);

struct vectorfunction
{
    int inputdim;  //输入维度
    int outputdim; //输出维度
    void (*function)(Vector *, Vector *);
};
typedef struct vectorfunction NdVectorfunction;
NdVectorfunction *ndVectorfunction_alloc(void (*function)(Vector *, Vector *), int intputdim, int outputdim);
void ndVectorfunction_free(NdVectorfunction *function);
void ndVectorfunction_call(NdVectorfunction *function, Vector *input, Vector *output);
void ndVectorfunction_jacobian(const NdVectorfunction *function, const Vector *x0, double h, Matrix *jacobian);