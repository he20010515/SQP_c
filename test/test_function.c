#include "function.h"
#include "stdio.h"
#include "vector.h"
double target_function(Vector *mat)
{
    return pow(mat->entry[0], 2) +
           pow(mat->entry[1], 2) +
           pow(mat->entry[2], 2);
}

int main(int argc, char const *argv[])
{
    printf("hello test function");
    NdsclaFunction *f = NdsclaFunctionAlloc(target_function, 3);

    return 0;
}