#include "lp.h"
#include "matrix.h"
#include "vector.h"
#include "util.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    Matrix *A = matrix_alloc(3, 5);
    double A_array[3][5] = {{
                                1,
                                0,
                                -1,
                                0,
                                0,
                            },
                            {
                                0,
                                1,
                                0,
                                -1,
                                0,
                            },
                            {
                                1,
                                1,
                                0,
                                0,
                                1,
                            }};
    array_2_matrix((double *)A_array, 3, 5, A);
    Vector *b = vector_alloc(3);
    b->entry[0] = 0;
    b->entry[1] = 0;
    b->entry[2] = 1;
    Vector *c = vector_alloc(5);
    c->entry[0] = 1;
    c->entry[1] = 1;
    c->entry[2] = 0;
    c->entry[3] = 0;
    c->entry[4] = 0;
    int init_base[3] = {2, 3, 4};
    Vector *x0 = vector_alloc(5);
    Vector *xstat = vector_alloc(5);
    optimize_lp_standard_type(c, b, A, x0, init_base, xstat);
    vector_print(xstat);
    return 0;
}
