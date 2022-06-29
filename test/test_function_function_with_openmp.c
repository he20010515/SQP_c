/*
 * @Author: HeYuwei
 * @Date: 2022-05-27 10:54:07
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-05-27 19:41:50
 * @FilePath: \SQP_c\test\test_function_function_with_openmp.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "function.h"
#include "matrix.h"
#include "math.h"
#include "function.h"
#include "vector.h"

double target_function(Vector *mat)
{
    float f = 0.0;
    for (size_t i = 0; i < mat->size; i++)
    {
        f += mat->entry[i] * mat->entry[i];
    }
    return f;
}

int main(void)
{
    sqp_init();
    NdsclaFunction *f = ndscla_function_alloc(target_function, 100);
    Vector *x0 = vector_alloc(100);
    vector_fill_const(x0, 100.0);
    Vector *grad = vector_alloc(100);
    ndscla_forward_grad(f, 0.01, x0, grad);
    // Matrix *hession = matrix_alloc(f->inputSize, f->inputSize);
    // ndscla_central_hession(f, 0.01, x0, hession);
    // matrix_print(hession);
}
