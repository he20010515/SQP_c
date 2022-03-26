#include "qp.h"
#include "matrix.h"
#include "vector.h"
#include "index_set.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    // Problme
    // f = (x-1)^2 + (y-2)^2 + (z-3)^2
    // s.t. x   == y
    //      x+z >= 2

    Matrix *G = matrix_alloc(3, 3);
    double G_array[3][3] = {
        {
            2,
            0,
            0,
        },
        {0, 2, 0},
        {0, 0, 2}};
    array_2_matrix((double *)G_array, 3, 3, G);

    Vector *c = vector_alloc(3);
    Vector *b = vector_alloc(2);
    c->entry[0] = -2.;
    c->entry[1] = -4;
    c->entry[2] = -6;

    b->entry[0] = 0.;
    b->entry[1] = 2.;

    Matrix *A = matrix_alloc(2, 3);
    double A_array[2][3] = {
        {1, -1, 0},
        {1, 0, 1},
    };
    array_2_matrix((double *)A_array, 2, 3, A);
    LinearConstraints *con = constraints_alloc(3, 2, 1, 1, A, b);
    Vector *x0 = vector_alloc(3);
    Vector *xstar = vector_alloc(3);
    x0->entry[0] = 6;
    x0->entry[1] = 6;
    x0->entry[2] = 6;

    optimize_qp_active_set(G, c, con, x0, xstar);
    vector_print(xstar);
    return 0;
}
