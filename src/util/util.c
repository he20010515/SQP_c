/*
 * @Author: HeYuwei
 * @Date: 2022-03-13 09:13:40
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-01 23:46:32
 * @FilePath: \SQP_c\src\util\util.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include <stdio.h>
#include <math.h>
#include "elog.h"
#include "util.h"
#include "vector.h"
#include "assert.h"
#include "matrix.h"
#include <stdlib.h>
#include <omp.h>

#define LOG_TAG "util"

void terminate(char *string)
{
    log_a("%s\n", string);
    log_a("The program is exiting now. . . .");
    fprintf(stdout, "The program is exiting now. . . .\n\n");
    Vector *v = vector_alloc(4);
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

void opem_mp_test(void)
{

    omp_set_num_threads(4);
#pragma omp parallel
    {
        printf("I'm thread %d\n", omp_get_thread_num());
    }
}