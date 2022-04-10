#include "simplex.h"
#include "util.h"
int main(int argc, char const *argv[])
{
    sqp_init();
    Vector *c = vector_alloc(3);
    c->entry[0] = 0.4;
    c->entry[1] = 0.5;
    c->entry[2] = 0;
    Matrix *A_ub = matrix_alloc(1, 3);
    A_ub->matrix_entry[0][0] = 0.3;
    A_ub->matrix_entry[0][1] = 0.1;
    A_ub->matrix_entry[0][2] = 0.0;
    Vector *b_ub = vector_alloc(1);
    b_ub->entry[0] = 2.7;
    Matrix *A_eq = matrix_alloc(2, 3);
    A_eq->matrix_entry[0][0] = 0.5;
    A_eq->matrix_entry[0][1] = 0.5;
    A_eq->matrix_entry[0][2] = 0.0;
    A_eq->matrix_entry[1][0] = 0.6;
    A_eq->matrix_entry[1][1] = 0.4;
    A_eq->matrix_entry[1][2] = -1.0;
    Vector *b_eq = vector_alloc(2);
    b_eq->entry[0] = 6.0;
    b_eq->entry[1] = 6.0;
    Vector *xstar = vector_alloc(3);
    _linprog_simplex(c, A_eq, b_eq, 100, 1e-8, FALSE, xstar);

    vector_print(xstar);

    return 0;
}
