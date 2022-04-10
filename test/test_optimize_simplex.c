#include "simplex.h"
#include "vector.h"
#include "matrix.h"
#include "lp.h"

int main(int argc, char const *argv[])
{
    // maxz = 3x1 -x2 -x3
    // s.t.
    //   x1 -2x2  +x3 <=11
    // -4x1  +x2 +2x3 >=3
    // -2x1       +x3 ==1
    double A_[3][3] = {
        {-2, 0, 1},
        {-4, 1, 2},
        {-1, 2, -1}};
    sqp_init();
    for (int i = 0; i < 100; i++)
    {
        Vector *b = vector_alloc(3);
        Matrix *A = matrix_alloc(3, 3);
        Vector *c = vector_alloc(3);
        LinearConstraints *con = linearconstraints_alloc(3, 3, 1, 2, A, b);
        Vector *x0 = vector_alloc(3);
        array_2_matrix((double *)A_, 3, 3, A);
        b->entry[0] = 1;
        b->entry[1] = 11;
        b->entry[2] = -3;
        c->entry[0] = -3;
        c->entry[1] = 1;
        c->entry[2] = 1;
        optimize_lp(con, c, x0);

        vector_print(x0);

        matrix_free(A);
        vector_free(c);
        vector_free(b);
        vector_free(x0);
        linearconstraints_free(con, 0);
    }
    return 0;
}
