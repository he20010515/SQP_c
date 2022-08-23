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
    printf("function->real_call_num:%ld\nfunction->recored_num:%ld\n", fun->real_call_num, fun->record_call_num);

    return 0;
}
