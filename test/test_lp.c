#include "lp.h"
#include "matrix.h"
#include "vector.h"

int main(int argc, char const *argv[])
{
    Matrix *A = matrix_alloc(3, 6);
    double A_array[3][6] = {{
                                2,
                                1,
                                0,
                                1,
                                0,
                                0,
                            },
                            {-4, -2, 3, 0, 1, 0},
                            {
                                1,
                                -2,
                                1,
                                0,
                                0,
                                1,
                            }};
    array_2_matrix((double *)A_array, 3, 6, A);
    Vector *b = vector_alloc(3);
    b->entry[0] = 8;
    b->entry[1] = 14;
    b->entry[2] = 18;
    Vector *c = vector_alloc(6);
    c->entry[0] = 6;
    c->entry[1] = -3;
    c->entry[2] = 1;
    c->entry[3] = 0;
    c->entry[4] = 0;
    c->entry[5] = 0;
    int init_base[3] = {3, 4, 5};
    Vector *x0 = vector_alloc(6);
    Vector *xstat = vector_alloc(6);
    optimize_lp_standard_type(c, b, A, init_base, xstat);
    vector_print(xstat);
    return 0;
}
