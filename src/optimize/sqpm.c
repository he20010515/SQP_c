/*
 * @Author: HeYuwei
 * @Date: 2022-07-15 12:48:08
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-18 15:14:11
 * @FilePath: \SQP_c\src\optimize\sqpm.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "matrix.h"
#include "vector.h"
#include "qp.h"
#include "function.h"
#include "sqp.h"
#include "elog.h"
#include <math.h>

void optimize_sqpm(const NdsclaFunction *fun,
                   const Nonlinearconstraints *con,
                   const Vector *x0,
                   const Vector *lambda0,
                   Vector *xstar)
{
    return;
}