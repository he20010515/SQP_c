#include <stdlib.h>
#include <stdio.h>
#include "vector.h"
#include <math.h>
#include "util.h"
#include "elog.h"

#define LOG_TAG "vector"
Vector *vector_alloc(int size)
{
    Vector *v = (Vector *)malloc(sizeof(Vector));
    v->entry = (double *)malloc(sizeof(double) * size);
    v->size = size;
    vector_fill_const(v, NAN);
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

void *vector_add_vector(const Vector *v, const Vector *w, Vector *v_w)
{
    if (v->size != w->size)
    {
        printf("v and w must have the same size");
        exit(-1);
    }

    for (size_t i = 0; i < v->size; i++)
    {
        v_w->entry[i] = v->entry[i] + w->entry[i];
    }
}

Vector *vector_multiply_const(const Vector *v, double a, int copy)
{
    if (copy == 1)
    {
        Vector *w = vector_alloc(v->size);
        for (size_t i = 0; i < v->size; i++)
        {
            w->entry[i] = v->entry[i] * a;
        }
        return w;
    }
    else
    {
        for (size_t i = 0; i < v->size; i++)
        {
            v->entry[i] *= a;
        }
        // return v;
        terminate("must copy");
    }
}

void *vector_print(const Vector *v)
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

int vector_argmin(const Vector *v)
{
    double min = NAN;
    int flag = 0;
    int iswaitfirstnumber = 1;
    for (int i = 0; i < v->size; i++)
    {
        if (!(isnan(v->entry[i]))AND iswaitfirstnumber) // 遇到的第一个非NAN的数设置0
        {
            min = v->entry[i];
            flag = i;
            iswaitfirstnumber = 0;
        }

        if (isnan(v->entry[i]))
            continue;

        if (v->entry[i] < min)
        {
            flag = i;
            min = v->entry[i];
        }
    }
    return flag;
}

int vector_argmax(const Vector *v)
{
    double min = NAN;
    int flag = 0;
    int iswaitfirstnumber = 1;
    for (int i = 0; i < v->size; i++)
    {
        if (!(isnan(v->entry[i]))AND iswaitfirstnumber) // 遇到的第一个非NAN的数设置0
        {
            min = v->entry[i];
            flag = i;
            iswaitfirstnumber = 0;
        }

        if (isnan(v->entry[i]))
            continue;

        if (v->entry[i] > min)
        {
            flag = i;
            min = v->entry[i];
        }
    }
    return flag;
}

double vector_max(const Vector *v)
{
    double min = NAN;
    int flag = 0;
    int iswaitfirstnumber = 1;
    for (int i = 0; i < v->size; i++)
    {
        if (!(isnan(v->entry[i]))AND iswaitfirstnumber) // 遇到的第一个非NAN的数设置0
        {
            min = v->entry[i];
            flag = i;
            iswaitfirstnumber = 0;
        }

        if (isnan(v->entry[i]))
            continue;

        if (v->entry[i] > min)
        {
            flag = i;
            min = v->entry[i];
        }
    }
    return min;
}
double vector_inner_product(const Vector *u, const Vector *v)
{
    if (!(u->size == v->size))
    {
        terminate("ERROR vector_inner_product size is not fit");
    }
    double sum = 0.0;
    for (int i = 0; i < u->size; i++)
    {
        sum = sum + u->entry[i] * v->entry[i];
    }
    return sum;
}

double vector_min(const Vector *v)
{
    double min = NAN;
    int flag = 0;
    int iswaitfirstnumber = 1;
    for (int i = 0; i < v->size; i++)
    {
        if (!isnan(v->entry[i]) AND iswaitfirstnumber) // 遇到的第一个非NAN的数设置0
        {
            min = v->entry[i];
            iswaitfirstnumber = 0;
        }

        if (isnan(v->entry[i]))
            continue;

        if (v->entry[i] < min)
        {
            flag = i;
            min = v->entry[i];
        }
    }
    return min;
}

double vector_1norm(const Vector *v)
{
    double m = 0;
    for (int i = 0; i < v->size; i++)
    {
        m = m + fabs(v->entry[i]);
    }
    return m;
}

void vector_fillna(Vector *v)
{
    for (int i = 0; i < v->size; i++)
    {
        if (isnan(v->entry[i]))
        {
            v->entry[i] = 0.0;
        }
    }
}
