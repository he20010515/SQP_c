/*
 * @Author: HeYuwei
 * @Date: 2022-05-20 18:36:52
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-05-27 10:44:03
 * @FilePath: \SQP_c\test\test_openmp.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include <stdio.h>

#include "omp.h"

int main(int argc, char const *argv[])
{
    int nthreads, tid;
    int share_int = 0;
#pragma omp parallel for num_threads(15) default(none) shared(share_int) private(nthreads, tid)
    for (int i = 0; i <= 10; i++)
    {
        tid = omp_get_thread_num();
        printf("hello world from thread = %d\n", tid);
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("number of threads = %d\n", nthreads);
        }
    }
    return 0;
}
