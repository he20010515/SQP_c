#include "function.h"
#include "matrix.h"
#include "math.h"
#include "function.h"
#include "vector.h"

#define M_alloc_variable(x, type, value)    \
    type *x = (type *)malloc(sizeof(type)); \
    *x = value

double target_function(Vector *mat)
{
    return pow(mat->entry[0], 2) +
           pow(mat->entry[1], 2) +
           pow(mat->entry[2], 2);
}

int main(void)
{
    NdsclaFunction *f = NdsclaFunctionAlloc(target_function, 3);
    Vector *x0 = vector_alloc(3);
    x0->entry[0] = 1.0;
    x0->entry[1] = 2.0;
    x0->entry[2] = 2.0;
    Vector *grad = vector_alloc(3);
    centralGrad(f, 0.01, x0, grad);
    vector_print(grad);
    printf("%f", NdsclaFunctionCall(f, x0));
}
