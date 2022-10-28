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
                   const Vector *mu0,
                   Vector *xstar)
{

    // // maxk=100;   %最大迭代次数
    // // n=length(x0); l=length(mu0); m=length(lam0);
    // // rho=0.5; eta=0.1;  B0=eye(n);
    // // x=x0; mu=mu0;  lam=lam0;
    // // Bk=B0; sigma=0.8;
    // // epsilon1=1e-6; epsilon2=1e-5;
    // // [hk,gk]=cons(x);  dfk=df1(x);
    // // [Ae,Ai]=dcons(x); Ak=[Ae; Ai];
    // // k=0;

    // size_t maxk = 100; // 最大迭代次数
    // int n = x0->size;
    // int l = con->i; // 不等式约束数量
    // int m = con->e; // 等式约束数量
    // double rho = 0.5;
    // double eta = 0.1;
    // double sigma = 0.8;
    // double epsilon1 = 1e-6;
    // double epsilon2 = 1e-5;

    // Matrix *B0 = matrix_callalloc(n);
    // Matrix *Bk = matrix_callalloc(n);
    // matrix_copy(B0, Bk);

    // Vector *x = vector_alloc(x->size);
    // vector_copy(x0, x);

    // Vector *mu = vector_alloc(mu0->size);
    // vector_copy(mu0, mu);

    // Vector *lambda = vector_alloc(lambda0->size);
    // vector_copy(lambda0, lambda);

    // int k = 0;

    // while (k <= maxk)
    // {
    //     ;
    // }

    // return;
}