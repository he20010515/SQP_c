/*
 * @Author: HeYuwei
 * @Date: 2022-04-09 18:42:11
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 18:15:53
 * @FilePath: \SQP_c\include\simplex.h
 * @Description: 仿照scipy中单纯形法的c语言版本
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#pragma once
#include "vector.h"
#include "matrix.h"
/**
 * @description:
 // Minimize a linear objective function subject to linear equality and
    // non-negativity constraints using the two phase simplex method.
    // Linear programming is intended to solve problems of the following form:

    // Minimize::

    //     c @ x

    // Subject to::

    //     A @ x == b
    //         x >= 0

    // User-facing documentation is in _linprog_doc.py.

    // Parameters
    // ----------
    // c : 1-D array
    //     Coefficients of the linear objective function to be minimized.
    // c0 : float
    //     Constant term in objective function due to fixed (and eliminated)
    //     variables. (Purely for display.)
    // A : 2-D array
    //     2-D array such that ``A @ x``, gives the values of the equality
    //     constraints at ``x``.
    // b : 1-D array
    //     1-D array of values representing the right hand side of each equality
    //     constraint (row) in ``A``.
    // Options
    // -------
    // maxiter : int
    //    The maximum number of iterations to perform.
    // tol : float
    //     The tolerance which determines when a solution is "close enough" to
    //     zero in Phase 1 to be considered a basic feasible solution or close
    //     enough to positive to serve as an optimal solution.
    // bland : bool
    //     If True, use Bland's anti-cycling rule [3]_ to choose pivots to
    //     prevent cycling. If False, choose pivots which should lead to a
    //     converged solution more quickly. The latter method is subject to
    //     cycling (non-convergence) in rare instances.

 * @param {Vector} *c
 * @param {Matrix} *A
 * @param {Vector} *b
 * @param {int} maxiter
 * @param {double} tol
 * @param {int} bland
 * @param {Vector} *x
 * @return {*}
 */
int _linprog_simplex(const Vector *c, const Matrix *A, const Vector *b, int maxiter, double tol, int bland, Vector *x);
