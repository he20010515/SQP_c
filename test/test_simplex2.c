#include "simplex.h"
#include "util.h"
int main(int argc, char const *argv[])
{
    sqp_init();
    //
    Vector *c = vector_alloc(4);
    c->entry[0] = -1;
    c->entry[1] = -1;
    c->entry[2] = 0;
    c->entry[3] = 0;
    Matrix *A_eq = matrix_alloc(2, 4);
    A_eq->matrix_entry[0][0] = 1.0;
    A_eq->matrix_entry[0][1] = 2.0;
    A_eq->matrix_entry[0][2] = 1.0;
    A_eq->matrix_entry[0][3] = 0.0;

    A_eq->matrix_entry[1][0] = 2.0;
    A_eq->matrix_entry[1][1] = 1.0;
    A_eq->matrix_entry[1][2] = 0.0;
    A_eq->matrix_entry[1][3] = 1.0;

    Vector *b_eq = vector_alloc(2);
    b_eq->entry[0] = 1.0;
    b_eq->entry[1] = 1.0;
    Vector *xstar = vector_alloc(4);
    Simplex *problem = simplex_alloc(c, NULL, NULL, A_eq, b_eq);
    simplex_main(problem, xstar);
    vector_print(xstar);

    return 0;
}
