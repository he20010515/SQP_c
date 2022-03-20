#pragma once
struct vector
{
    int size;
    double *entry;
};
typedef struct vector Vector;

Vector *vector_alloc(int size);
Vector *vector_free(Vector *v);
void vector_copy(const Vector *v, Vector *w);
Vector *vector_add_const(Vector *v, double a, int copy);
Vector *vector_add_vector(const Vector *v, const Vector *w, Vector *v_w);
Vector *vector_multiply_const(Vector *v, double a, int copy);
void *vector_print(Vector *v);
double vector_2norm(const Vector *v);
double vector_2metric(const Vector *v, const Vector *u);
void vector_fill_const(Vector *v, double a);
int vector_argmin(const Vector *v);
double vector_inner_product(const Vector *u, const Vector *v);
double vector_min(const Vector *v);