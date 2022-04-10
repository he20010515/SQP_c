#include "function.h"
#include "matrix.h"
#include "math.h"
#include "function.h"
#include "vector.h"

double target_function(Vector *mat)
{
    return pow(mat->entry[1], 2) * pow(mat->entry[2], 2) +
           pow(mat->entry[0], 2) * pow(mat->entry[2], 2) +
           pow(mat->entry[0], 2) * pow(mat->entry[1], 2);
}

int main(void)
{
    sqp_init();
    NdsclaFunction *f = ndscla_function_alloc(target_function, 3);
    Vector *x0 = vector_alloc(3);
    x0->entry[0] = 1.0;
    x0->entry[1] = 2.0;
    x0->entry[2] = 3.0;
    Vector *grad = vector_alloc(3);
    ndscla_central_grad(f, 0.01, x0, grad);
    vector_print(grad);
    Matrix *hession = matrix_alloc(f->inputSize, f->inputSize);
    ndscla_central_hession(f, 0.01, x0, hession);
    matrix_print(hession);
}
