#include "simplex.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    double A_array[5][9] = {
        {-1, 2, 1, 2, 1, 0, 0, 0, 0},
        {1, 2, -1, -2, 0, 1, 0, 0, 0},
        {1, -2, -1, 2, 0, 0, 1, 0, 0},
        {1, 0, -1, 0, 0, 0, 0, -1, 0},
        {0, 1, 0, -1, 0, 0, 0, 0, -1},
    };
    double beq[5] = {2, 6, 2, 0, 0};
    double c[9] = {-2, -5, 2, 5, 0, 0, 0, 0, 0};
    Vector b;
    Vector cc;
    Matrix *mat = matrix_alloc(5, 9);
    array_2_matrix((double *)A_array, 5, 9, mat);

    b.entry = beq;
    b.size = 5;

    cc.entry = c;
    cc.size = 9;
    Vector *x = vector_alloc(9);
    _linprog_simplex(&cc, mat, &b, 30, 1e-5, FALSE, x);
    vector_print(x);
    return 0;
}
