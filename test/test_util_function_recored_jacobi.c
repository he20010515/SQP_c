/*
 * @Author: heyuwei he20010515@163.com
 * @Date: 2022-09-18 18:59:52
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-13 23:42:54
 * @FilePath: /SQP_c/test/test_util_function_recored.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "function.h"
#include "math.h"
void test_function(const Vector *x, Vector *y)
{
    for (int j = 0; j < y->size; j++)
    {
        double p = 0;
        for (int i = 0; i < x->size; i++)
        {
            if (i == j)
            {
                break;
            }
            else
            {
                p += x->entry[i] * x->entry[i];
            }
        }
        y->entry[j] = p;
    }
}

int main(int argc, char const *argv[])
{
    const int N = 10;
    const int M = 10;
    NdVectorfunction *fun = ndVectorfunction_alloc(test_function, N, M);
    Vector *x0 = vector_alloc(N);
    for (int i = 0; i < N; i++)
    {
        x0->entry[i] = 2;
    }
    Matrix *jacobian = matrix_alloc(M, N);
    ndVectorfunction_jacobian(fun, x0, 0.2, jacobian);
    matrix_print(jacobian);

    ndVectorfunction_free(fun);
    vector_free(x0);

    return 0;
}
