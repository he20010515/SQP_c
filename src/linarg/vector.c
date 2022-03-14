#include <stdlib.h>
#include <stdio.h>
#include "vector.h"
#include <math.h>
#include "util.h"
Vector *vector_alloc(int size)
{
    Vector *v = (Vector *)malloc(sizeof(Vector));
    v->entry = (double *)malloc(sizeof(double) * size);
    v->size = size;
    return v;
};

Vector *vector_free(Vector *v)
{
    free(v->entry);
    free(v);
}

void vector_copy(const Vector *v, Vector *w)
{
    // copy v to w
    if (v->size != w->size)
    {
        printf("vector_copy: w and v must have the same size");
        exit(-1);
    }

    for (size_t i = 0; i < v->size; i++)
    {
        w->entry[i] = v->entry[i];
    }
    w->size = v->size;
    return;
}

Vector *vector_add_const(Vector *v, double a, int copy)
{
    if (copy == 1)
    {
        Vector *w = vector_alloc(v->size);
        vector_copy(v, w);
        for (size_t i = 0; i < v->size; i++)
        {
            w->entry[i] += a;
        }
        return w;
    }

    for (size_t i = 0; i < v->size; i++)
    {
        v->entry[i] += a;
    }
    return v;
}

Vector *vector_add_vector(const Vector *v, const Vector *w, Vector *v_w)
{
    if (v->size != w->size)
    {
        printf("v and w must have the same size");
        exit(-1);
    }

    Vector *u = vector_alloc(v->size);
    for (size_t i = 0; i < u->size; i++)
    {
        v_w->entry[i] = v->entry[i] + w->entry[i];
    }
}

Vector *vector_multiply_const(Vector *v, double a, int copy)
{
    if (copy == 1)
    {
        Vector *w = vector_alloc(v->size);
        for (size_t i = 0; i < v->size; i++)
        {
            w->entry[i] = v->entry[i] * a;
        }
    }
    else
    {
        for (size_t i = 0; i < v->size; i++)
        {
            v->entry[i] *= a;
        }
    }
}

void *vector_print(Vector *v)
{
    printf("[");
    for (size_t i = 0; i < v->size; i++)
    {
        printf("%g, ", v->entry[i]);
    }
    printf("]\n");
}

double vector_2norm(const Vector *v)
{
    double s = 0.0;
    for (size_t i = 0; i < v->size; i++)
    {
        s += pow(v->entry[i], 2);
    }
    return sqrt(s);
}

double vector_2metric(const Vector *v, const Vector *u)
{
    if (v->size != u->size)
    {
        terminate("ERROR: vector_2metric: vector must have the same size");
    }

    Vector *temp = vector_alloc(v->size);
    for (size_t i = 0; i < v->size; i++)
    {
        temp->entry[i] = v->entry[i] - u->entry[i];
    }
    double distance = vector_2norm(temp);
    vector_free(temp);
    return distance;
}

void vector_fill_const(Vector *v, double a)
{
    for (size_t i = 0; i < v->size; i++)
    {
        v->entry[i] = a;
    }
    return;
}