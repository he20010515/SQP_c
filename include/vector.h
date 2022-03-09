#pragma once

struct vector
{
    int size;
    double *entry;
};
typedef struct vector Vector;

Vector *vector_alloc(int size);
Vector *vector_free(Vector *v);
void vector_copy(Vector *v, Vector *w);
Vector *vector_add_const(Vector *v, double a, int copy);
Vector *vector_add_vector(Vector *v, Vector *w, int copy);
Vector *vector_multiply_const(Vector *v, double a, int copy);
void *vector_print(Vector *v);