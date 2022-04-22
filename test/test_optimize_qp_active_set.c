#include "qp.h"
#include "matrix.h"
#include "vector.h"
#include "index_set.h"
#include "util.h"
int main(int argc, char const *argv[])
{
    sqp_init();

    double G_array[2][2] = {2, 0, 0, 2};
    double A_array[5][2] = {
        {
            1,
            -2,
        },
        {
            -1,
            -2,
        },
        {-1, 2},
        {1, 0},
        {0, 1},
    };
    for (int i = 0; i < 5; i++)
    {
        Matrix *G = matrix_alloc(2, 2);
        array_2_matrix((double *)G_array, 2, 2, G);
        Vector *c = vector_alloc(2);
        Vector *b = vector_alloc(5);

        c->entry[0] = -2.;
        c->entry[1] = -5.;
        b->entry[0] = -2.;
        b->entry[1] = -6.;
        b->entry[2] = -2.;
        b->entry[3] = 0.;
        b->entry[4] = 0.;

        Matrix *A = matrix_alloc(5, 2);
        array_2_matrix((double *)A_array, 5, 2, A);
        LinearConstraints *con = linearconstraints_alloc(2, 5, 0, 5, A, b);
        Vector *x0 = vector_alloc(2);
        Vector *xstar = vector_alloc(2);
        x0->entry[0] = 2;
        x0->entry[1] = 0;
        Vector *lambda = vector_alloc(5);
        optimize_qp_active_set(G, c, con, NULL, xstar, lambda);
        vector_print(xstar);
        vector_print(lambda);

        matrix_free(G);
        vector_free(c);
        vector_free(b);

        matrix_free(A);
        vector_free(x0);
        vector_free(xstar);
        vector_free(lambda);
        linearconstraints_free(con, FALSE);
    }

    return 0;
}
