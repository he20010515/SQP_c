#include "linear_equations.h"
#include "matrix.h"
#include "vector.h"

int main(int argc, char const *argv[])
{
    double A_aray[4][4] = {
        {1, 0, -1, 0},
        {0, 1, 2, 1},
        {-1, 2, 0, 0},
        {0, 1, 0, 0},
    };
    Vector *b = vector_alloc(4);
    b->entry[0] = 1;
    b->entry[1] = 2;
    b->entry[2] = 3;
    b->entry[3] = 4;

    Matrix *A = matrix_alloc(4, 4);
    array_2_matrix((double *)A_aray, 4, 4, A);
    Vector *x = vector_alloc(4);
    linear_equation_gaussian_elimination(A, b, x);
    vector_print(x);
    return 0;
}
