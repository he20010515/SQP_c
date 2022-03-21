#include <stdio.h>
#include <math.h>
#include "elog.h"
#include "util.h"
#include "vector.h"
#include "assert.h"
#include <stdlib.h>

#define LOG_TAG "util"
void terminate(char *string)
{
    log_a("%s\n", string);
    log_a("The program is exiting now. . . .");
    fprintf(stdout, "The program is exiting now. . . .\n\n");
    Vector *v = vector_alloc(4);
    v->entry[4200] = 42.0;
    exit(-1);
}

int double_equal(double a, double b)
{
    if (fabs(a - b) < DOUBLE_ERROR)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int vector_any_bigger_equal_than_const(const Vector *v, double a)
{
    for (int i = 0; i < v->size; i++)
    {
        if (v->entry[i] < a)
            return 0;
    }
    return 1;
}