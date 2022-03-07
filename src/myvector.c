#include <stdlib.h>
#include <stdio.h>
#include "myvector.h"

Vector *VectorAlloc(int size)
{
    Vector *v = (Vector *)malloc(sizeof(Vector));
    v->entry = (double *)malloc(sizeof(double) * size);
    v->size = size;
};

Vector *VectorFree(Vector *v)
{
    free(v->entry);
    free(v);
}

void VectorCopy(Vector *v, Vector *w)
{
    if (v->size != w->size)
    {
        printf("VectorCopy: w and v must have the same size");
        exit(-1);
    }

    for (size_t i = 0; i < v->size; i++)
    {
        w->entry[i] = v->entry[i];
    }
    w->size = v->size;
}

Vector *VectorAddConst(Vector *v, double a, int copy)
{
    if (copy == 1)
    {
        Vector *w = VectorAlloc(v->size);
        VectorCopy(v, w);
        for (size_t i = 0; i < v->size; i++)
        {
            w->entry[i] += a;
        }
    }

    for (size_t i = 0; i < v->size; i++)
    {
        v->entry[i] += a;
    }
}

Vector *VectorAddVector(Vector *v, Vector *w, int copy)
{
    if (v->size != w->size)
    {
        printf("v and w must have the same size");
        exit(-1);
    }
    if (copy == 1)
    {
        Vector *u = VectorAlloc(v->size);
        for (size_t i = 0; i < u->size; i++)
        {
            u->entry[i] = v->entry[i] + w->entry[i];
        }
        return u;
    }
    else
    {
        for (size_t i = 0; i < v->size; i++)
        {
            w->entry[i] += v->entry[i];
        }
        return w;
    }
}

Vector *VectorMultiply(Vector *v, double a, int copy)
{
    if (copy == 1)
    {
        Vector *w = VectorAlloc(v->size);
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

void *VectorPrint(Vector *v)
{
    printf("[");
    for (size_t i = 0; i < v->size; i++)
    {
        printf("%f, ", v->entry[i]);
    }
    printf("]");
}