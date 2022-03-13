#include "matrix.h"
Matrix *foo()
{
    Matrix *mat = matrix_alloc(5, 1);
    return mat;
}

int main(int argc, char const *argv[])
{
    Matrix *mat = matrix_alloc(5, 1);
    matrix_free(mat);

    Matrix *matrix_from_foo = foo();
    matrix_free(matrix_from_foo);
    return 0;
}