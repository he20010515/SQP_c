/*
 * @Author: heyuwei he20010515@163.com
 * @Date: 2022-09-18 18:59:52
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-13 22:27:03
 * @FilePath: /SQP_c/test/test_util_function_recored.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "function.h"

double test_function(Vector *x)
{
    double y = 0.0;
    for (int i = 0; i < x->size; i++)
    {
        y += x->entry[i];
    }

    return y;
}

int main(int argc, char const *argv[])
{

    NdsclaFunction *fun = ndscla_function_alloc(test_function, 5);
    Vector *x0 = vector_alloc(5);
    for (int i = 0; i < 4; i++)
    {
        x0->entry[0] = i;
        x0->entry[1] = i + 1;
        x0->entry[2] = i + 2;
        x0->entry[3] = i + 3;
        x0->entry[4] = i + 4;
        double temp = 0.0;
        temp = ndscla_function_call(fun, x0);
        printf("x[%d],f(x[i]) = %lf\n", i, temp);
    }

    for (int i = 0; i < 4; i++)
    {
        x0->entry[0] = i;
        x0->entry[1] = i + 1;
        x0->entry[2] = i + 2;
        x0->entry[3] = i + 3;
        x0->entry[4] = i + 4;
        double temp = 0.0;
        temp = ndscla_function_call(fun, x0);
        printf("times2 =x[%d],f(x[i]) = %lf\n", i, temp);
    }
    ndscla_function_free(fun);
    vector_free(x0);
    printf("function->real_call_num:%ld\nfunction->recored_num:%ld\n", fun->real_call_num, fun->record_call_num);

    return 0;
}
