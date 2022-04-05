#include "qp.h"
#include "matrix.h"
#include "vector.h"
#include "index_set.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    Matrix *G = matrix_alloc(2, 2);
    double G_array[2][2] = {2, -2, -2, 4};
    array_2_matrix((double *)G_array, 2, 2, G);

    Vector *c = vector_alloc(2);
    Vector *b = vector_alloc(5);
    c->entry[0] = -4.;
    c->entry[1] = -12.;
    b->entry[0] = 2.;
    b->entry[1] = -2.;
    b->entry[2] = -3.;
    b->entry[3] = 0.;
    b->entry[4] = 0.;

    Matrix *A = matrix_alloc(5, 2);
    double A_array[5][2] = {
        {
            1,
            1,
        },
        {
            1,
            -2,
        },
        {-2, -1},
        {1, 0},
        {0, 1},
    };
    array_2_matrix((double *)A_array, 5, 2, A);
    LinearConstraints *con = linearconstraints_alloc(2, 5, 1, 4, A, b);
    Vector *x0 = vector_alloc(2);
    Vector *lambda = vector_alloc(5);
    Vector *xstar = vector_alloc(2);
    x0->entry[0] = 1.2;
    x0->entry[1] = 0.8;
    optimize_qp_active_set(G, c, con, x0, xstar, lambda);
    vector_print(xstar);
    return 0;
}
