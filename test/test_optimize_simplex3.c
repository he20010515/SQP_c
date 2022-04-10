#include "simplex.h"
#include "vector.h"
#include "matrix.h"
#include "lp.h"
#include "random.h"

int main(int argc, char const *argv[])
{
    double A_[3][3] = {
        {-2, 0, 1},
        {-4, 1, 2},
        {-1, 2, -1}};
    const int x = 4000;
    const int y = 5 * x;
    sqp_init();
    for (int i = 0; i < 1; i++)
    {
        Vector *b = vector_alloc(x);
        Matrix *A = matrix_alloc(x, y);
        Vector *c = vector_alloc(y);

        LinearConstraints *con = linearconstraints_alloc(y, x, x, 0, A, b);
        Vector *x0 = vector_alloc(3);
        for (int i = 0; i < A->row_size; i++)
        {
            for (int j = 0; j < A->col_size; j++)
            {
                A->matrix_entry[i][j] = rand_int(0, 20);
            }
            c->entry[i] = rand_int(-10, 20);
        }
        for (int j = 0; j < b->size; j++)
        {
            b->entry[i] = rand_int(0, 50);
        }
        optimize_lp(con, c, x0,100,1e-7,FALSE);
        // vector_print(x0);
        matrix_free(A);
        vector_free(c);
        vector_free(b);
        vector_free(x0);
        linearconstraints_free(con, 0);
    }
    return 0;
}
