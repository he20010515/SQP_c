// /*
//  * @Author: HeYuwei
//  * @Date: 2022-07-15 12:48:08
//  * @LastEditors: Heyuwei
//  * @LastEditTime: 2022-07-18 15:14:11
//  * @FilePath: \SQP_c\src\optimize\sqpm.c
//  * @Description:
//  *
//  * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
//  */

#include "matrix.h"
// #include "vector.h"
// #include "qp.h"
// #include "function.h"
// #include "sqp.h"
// #include "elog.h"
// #include <math.h>

// void qsubp(Vector *vector_dk,
//            Vector *vector_mu,
//            Vector *vector_lam,
//            int kdimension_x,
//            int knum_eq_constraints,
//            int knum_ineq_constraints,
//            Vector *vector_Df,
//            int matrix_Bk,
//            Vector *vector_h_eq,
//            Vector *vector_g_ineq,
//            Matrix *matrix_Dh_eq,
//            Matrix *matrix_Dg_ineq)
// {
//     Vector *CON_RHS = vector_alloc(knum_eq_constraints + knum_ineq_constraints); // 约束的右边项
//     Matrix *CON_LHS = matrix_alloc(knum_eq_constraints + knum_ineq_constraints, kdimension_x);
//     Vector *lambda = vector_alloc(knum_eq_constraints + knum_ineq_constraints);

//     for (int i = 0; i < knum_eq_constraints; i++)
//         CON_RHS->entry[i] = vector_h_eq->entry[i];

//     for (int i = 0; i < knum_ineq_constraints; i++)
//         CON_RHS->entry[i + knum_eq_constraints] = vector_g_ineq->entry[i];

//     for (int j = 0; j < kdimension_x; j++)
//     {
//         for (int i = 0; i < knum_eq_constraints; i++)
//             CON_LHS->matrix_entry[i][j] = matrix_Dh_eq->matrix_entry[i][j];
//         for (int i = 0; i < knum_ineq_constraints; i++)
//             CON_LHS->matrix_entry[i + knum_ineq_constraints][j] = matrix_Dg_ineq->matrix_entry[i][j];
//     }

//     LinearConstraints *cons = linearconstraints_alloc(kdimension_x, knum_eq_constraints + knum_ineq_constraints, knum_eq_constraints, knum_ineq_constraints, CON_LHS, CON_RHS);
//     optimize_qp_active_set(matrix_Bk, vector_Df, cons, NULL, vector_dk, lambda);
//     for (int i = 0; i < knum_eq_constraints; i++)
//         vector_mu->entry[i] = lambda->entry[i];
//     for (int i = 0; i < knum_ineq_constraints; i++)
//         vector_lam->entry[i] = lambda->entry[i + knum_eq_constraints];

//     vector_free(CON_RHS);
//     matrix_free(CON_LHS);
//     vector_free(lambda);
//     linearconstraints_free(cons, FALSE);
//     return;
// }

// double max_(double a, double b)
// {
//     if (a > b)
//     {
//         return a;
//     }
//     else
//     {
//         return b;
//     }
// }

// double Value_Fun(double sigma, double obj_fun, double norm_eq_constraint, double norm_ineq_constraint, const int knum_eq_constraints, const int knum_ineq_constraints)
// {
//     double Value_Fun_result;
//     if (knum_eq_constraints == 0)
//     {
//         Value_Fun_result = obj_fun + 1.0 / sigma * norm_ineq_constraint;
//     }
//     else if (knum_ineq_constraints == 0)
//     {
//         Value_Fun_result = obj_fun + 1.0 / sigma * norm_eq_constraint;
//     }
//     else if (knum_eq_constraints > 0 && knum_ineq_constraints > 0)
//     {
//         Value_Fun_result = obj_fun + 1.0 / sigma * (norm_eq_constraint + norm_ineq_constraint);
//     }

//     return Value_Fun_result;
// }
// double Jacobi_Value_Fun(double sigama, Vector *matrix_d, Vector *vector_Df, double norm_eq_constraint, double norm_ineq_constraint, const int knum_eq_constraints, const int knum_ineq_constraints)
// {
//     // matrix_d should be column vectors
//     // printf("    Start calculating the Jacobi matrix of Lagrangian Function\n");

//     double multiply_result = vector_inner_product(vector_Df, matrix_d);

//     // double Jacobi_Value_Fun_result;
//     if (knum_eq_constraints == 0)
//     {
//         return multiply_result - 1.0 / sigama * norm_ineq_constraint;
//     }
//     else if (knum_ineq_constraints == 0)
//     {
//         return multiply_result - 1.0 / sigama * norm_eq_constraint;
//     }
//     else if (knum_eq_constraints > 0 && knum_ineq_constraints > 0)
//     {
//         return multiply_result - 1.0 / sigama * (norm_eq_constraint + norm_ineq_constraint);
//     }

//     // return Jacobi_Value_Fun_result;
// }

// void vector_a_add_k_b(const Vector *a, const double k, const Vector *b, Vector *res)
// {
//     for (int i = 0; i < a->size; i++)
//         res->entry[i] = a->entry[i] + k * b->entry[b];

//     return;
// }

// void optimize_sqpm(const NdsclaFunction *Objective_Fun,
//                    const Nonlinearconstraints *con,
//                    const Vector *x_init,
//                    const Vector *mu_init,
//                    const Vector *lam_init,
//                    Vector *xstar)
// {
//     const int kdimension_x = x_init->size;
//     const int knum_eq_constraints = con->e;
//     const int knum_ineq_constraints = con->i;

//     const int kmax_iteration = 200;

//     // Standard search :
//     //			 m : number of searches
//     //		m_stop : number of times the search was stopped
//     //		 max_m : maximum number of searches
//     //	   relax_m : 累计松弛搜索后的迭代次数
//     // max_relax_m : 松弛搜索后 5 次迭代目标函数没有下降则返回松弛搜索前的点

//     int i, num_iterations = 0, num_iteration_min = 0, num_iteration_relax = 0,
//            m = 0, m_stop = 0, max_m = 8, mid_m = 4, relax_m = 0, max_relax_m = 5;
//     double rho = 0.3,
//            gamma = 0.5,
//            eta = 0.1,
//            tau = 0,
//            sigma = 1,
//            epsilon1 = 1e-3,
//            epsilon2 = 1e-2,
//            epsilon3 = 1e-3,
//            deta = 0.05,
//            theta = 0;

//     //! 一个用于组合两个约束变量的临时向量:
//     Vector *vector_temp_c = vector_alloc(knum_eq_constraints + knum_ineq_constraints);
//     Matrix *matrix_temp_jacobi = matrix_alloc(kdimension_x, knum_eq_constraints + knum_ineq_constraints);

//     // vector_alloc():
//     Vector *vector_x = vector_alloc(kdimension_x);
//     Vector *vector_next_x = vector_alloc(kdimension_x);
//     Vector *vector_min_x = vector_alloc(kdimension_x);
//     Vector *vector_relax_x = vector_alloc(kdimension_x);
//     Vector *vector_mu = vector_alloc(knum_eq_constraints);
//     Vector *vector_min_mu = vector_alloc(knum_eq_constraints);
//     Vector *vector_relax_mu = vector_alloc(knum_eq_constraints);
//     Vector *vector_lam = vector_alloc(knum_ineq_constraints);
//     Vector *vector_min_lam = vector_alloc(knum_ineq_constraints);
//     Vector *vector_relax_lam = vector_alloc(knum_ineq_constraints);
//     Vector *vector_Df = vector_alloc(kdimension_x);
//     Vector *vector_next_Df = vector_alloc(kdimension_x);
//     Vector *vector_min_Df = vector_alloc(kdimension_x);
//     Vector *vector_relax_Df = vector_alloc(kdimension_x);
//     Vector *vector_h_eq = vector_alloc(knum_eq_constraints);
//     Vector *vector_next_h_eq = vector_alloc(knum_eq_constraints);
//     Vector *vector_min_h_eq = vector_alloc(knum_eq_constraints);
//     Vector *vector_relax_h_eq = vector_alloc(knum_eq_constraints);
//     Vector *vector_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_next_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_min_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_relax_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_neg_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_next_neg_g_ineq = vector_alloc(knum_ineq_constraints);
//     Vector *vector_mu_lam = vector_alloc(knum_eq_constraints + knum_ineq_constraints);
//     Vector *vector_dk = vector_alloc(kdimension_x);
//     Vector *vector_relax_dk = vector_alloc(kdimension_x);

//     Matrix *matrix_Bk = matrix_alloc(kdimension_x, kdimension_x);
//     Matrix *matrix_relax_Bk = matrix_alloc(kdimension_x, kdimension_x);
//     Matrix *matrix_Dh_Dg = matrix_alloc(knum_eq_constraints + knum_ineq_constraints, kdimension_x);
//     Matrix *matrix_Dh_eq = matrix_alloc(knum_eq_constraints, kdimension_x);
//     Matrix *matrix_next_Dh_eq = matrix_alloc(knum_eq_constraints, kdimension_x);
//     Matrix *matrix_min_Dh_eq = matrix_alloc(knum_eq_constraints, kdimension_x);
//     Matrix *matrix_relax_Dh_eq = matrix_alloc(knum_eq_constraints, kdimension_x);
//     Matrix *matrix_Dg_ineq = matrix_alloc(knum_ineq_constraints, kdimension_x);
//     Matrix *matrix_next_Dg_ineq = matrix_alloc(knum_ineq_constraints, kdimension_x);
//     Matrix *matrix_min_Dg_ineq = matrix_alloc(knum_ineq_constraints, kdimension_x);
//     Matrix *matrix_relax_Dg_ineq = matrix_alloc(knum_ineq_constraints, kdimension_x);

//     matrix_fill_const(matrix_Bk, 0.0);
//     for (int i = 0; i < matrix_Bk->row_size; i++)
//         matrix_Bk->matrix_entry[i][i] = 1.0;

//     double
//         obj_fun,
//         obj_fun_next, obj_fun_min, obj_fun_relax,
//         value_fun, value_fun_next, value_fun_min, value_fun_relax,
//         jacobi_value_fun, jacobi_value_fun_min, jacobi_value_fun_relax;

//     // search step
//     double alpha = rho;

//     Vector *penalty_parameter_1 = vector_alloc(kdimension_x);
//     Vector *penalty_parameter_2 = vector_alloc(kdimension_x);
//     Vector *penalty_parameter_min = vector_alloc(kdimension_x);
//     Vector *penalty_parameter_min_add1 = vector_alloc(kdimension_x);

//     // 用于 Armijo 搜索的中间变量
//     double differ_value_fun = 0, value_fun_1 = 0, value_fun_2 = 0;
//     // 1 标准
//     // 0 松弛

//     // flag of optimal stopping criterion 1
//     double norm_eq_constraint, norm_next_eq_constraint,
//         norm_ineq_constraint, norm_next_ineq_constraint;

//     // flag of optimal stopping criterion 2
//     bool Flag_break = TRUE;
//     bool flag_search = TRUE;
//     bool flag_relax = FALSE;
//     // flag_search: Choose to use standard or relaxed search
//     // Standard search for the first time

//     // initialize
//     vector_copy(x_init, vector_x);
//     vector_copy(x_init, vector_next_x);
//     vector_copy(x_init, vector_min_x);

//     vector_copy(mu_init, vector_mu);
//     vector_copy(mu_init, vector_min_mu);
//     vector_copy(lam_init, vector_lam);
//     vector_copy(lam_init, vector_min_lam);

//     obj_fun = ndscla_function_call(Objective_Fun, vector_x);

//     ndscla_forward_grad(Objective_Fun, 0, vector_x, vector_Df); // TIP 这里的差分步长是自动设定的,所以传进去0

//     ndVectorfunction_call(con->c, vector_x, vector_temp_c);
//     for (int i = 0; i < knum_eq_constraints; i++)
//         vector_h_eq->entry[i] = vector_temp_c->entry[i];
//     for (int i = 0; i < knum_ineq_constraints; i++)
//         vector_g_ineq->entry[i] = vector_temp_c->entry[i + knum_eq_constraints];

//     ndVectorfunction_jacobian(con->c, vector_x, h, matrix_temp_jacobi); // TIP 这里的差分步长需要手动设定
//     for (int i = 0; i < knum_eq_constraints; i++)
//         for (int j = 0; j < kdimension_x; j++)
//             matrix_Dh_eq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i][j];
//     for (int i = 0; i < knum_ineq_constraints; i++)
//         for (int j = 0; j < kdimension_x; j++)
//             matrix_Dg_ineq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i + knum_eq_constraints][j];

//     obj_fun_min = obj_fun;

//     vector_copy(vector_Df, vector_min_Df);
//     vector_copy(vector_h_eq, vector_min_h_eq);
//     vector_copy(vector_g_ineq, vector_min_g_ineq);

//     value_fun_min = INFINITY;
//     jacobi_value_fun_min = INFINITY;

//     while (Flag_break AND(num_iterations < kmax_iteration))
//     {
//         // calculating the subproblem
//         vector_fill_const(vector_mu, 0.0);
//         vector_fill_const(vector_lam, 0.0);
//         vector_fill_const(vector_dk, 1.0);
//         qsubp(vector_dk, vector_mu, vector_lam, kdimension_x, knum_eq_constraints, knum_ineq_constraints, vector_Df, matrix_Bk, vector_h_eq, vector_g_ineq, matrix_Dh_eq, matrix_Dg_ineq);
//         for (int i = 0; i < knum_ineq_constraints; i++)
//             vector_neg_g_ineq->entry[i] = max_(-1 * vector_g_ineq->entry[i], 0.0);

//         norm_eq_constraint = vector_1norm(vector_h_eq);
//         norm_ineq_constraint = vector_1norm(vector_neg_g_ineq);

//         value_fun = Value_Fun(sigma, obj_fun, norm_eq_constraint, norm_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);
//         jacobi_value_fun = Jacobi_Value_Fun(sigma, vector_dk, vector_Df, norm_eq_constraint, norm_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);

//         tau = max_(vector_max(vector_mu), vector_max(vector_lam)); // TODU debug here
//         if (sigma * (tau + deta) > 1.0)
//         {
//             sigma = 1.0 / (tau + 2.0 * deta);
//         }
//         if (flag_search)
//         {
//             // Standard search
//             m = 0;
//             m_stop = 0;
//             while (m < max_m)
//             {
//                 alpha = pow(rho, m);

//                 vector_a_add_k_b(vector_x, alpha, vector_dk, vector_next_x);

//                 obj_fun_next = ndscla_function_call(Objective_Fun, vector_next_x);

//                 ndVectorfunction_call(con->c, vector_x, vector_temp_c);
//                 for (int i = 0; i < knum_eq_constraints; i++)
//                     vector_next_h_eq->entry[i] = vector_temp_c->entry[i];
//                 for (int i = 0; i < knum_ineq_constraints; i++)
//                     vector_next_h_eq->entry[i] = vector_temp_c->entry[i + knum_eq_constraints];

//                 for (int i = 0; i < knum_ineq_constraints; i++)
//                     vector_neg_g_ineq->entry[i] = max_(-1 * vector_next_g_ineq->entry[i], 0);

//                 norm_next_eq_constraint = vector_1norm(vector_next_h_eq);
//                 norm_next_ineq_constraint = vector_1norm(vector_next_neg_g_ineq);

//                 value_fun_next = Value_Fun(sigma, obj_fun_next, norm_next_eq_constraint, norm_next_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);

//                 if (value_fun_next - value_fun < eta * alpha * jacobi_value_fun)
//                 {
//                     m_stop = m;
//                     break;
//                 }
//                 m += 1;
//             }
//             if (m == max_m)
//             {
//                 m_stop = mid_m;
//                 alpha = pow(rho, m_stop);
//                 vector_a_add_k_b(vector_x, alpha, vector_dk, vector_next_x);
//                 obj_fun_next = ndscla_function_call(Objective_Fun, vector_next_x);
//                 ndVectorfunction_call(con->c, vector_x, vector_temp_c);
//                 for (int i = 0; i < knum_eq_constraints; i++)
//                     vector_next_h_eq->entry[i] = vector_temp_c->entry[i];
//                 for (int i = 0; i < knum_ineq_constraints; i++)
//                     vector_next_g_ineq->entry[i] = vector_temp_c->entry[i + knum_eq_constraints];
//                 for (int i = 0; i < knum_ineq_constraints; i++)
//                     vector_neg_g_ineq->entry[i] = max_(-1 * vector_next_g_ineq->entry[i], 0);

//                 norm_next_eq_constraint = vector_1norm(vector_next_h_eq);
//                 norm_next_ineq_constraint = vector_1norm(vector_next_neg_g_ineq);

//                 value_fun_next = Value_Fun(sigma, obj_fun_next, norm_next_eq_constraint, norm_next_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);
//             }
//         }
//         else
//         {
//             // Relaxes search
//             {
//                 num_iteration_relax = num_iterations;
//                 obj_fun_relax = obj_fun;
//                 vector_copy(vector_x, vector_relax_x);
//                 vector_copy(vector_Df, vector_relax_Df);
//                 vector_copy(vector_h_eq, vector_relax_h_eq);
//                 vector_copy(vector_g_ineq, vector_relax_g_ineq);

//                 matrix_copy(matrix_Dh_eq, matrix_relax_Dh_eq);
//                 matrix_copy(matrix_Dg_ineq, matrix_relax_Dg_ineq);
//                 matrix_copy(matrix_Bk, matrix_relax_Bk);

//                 vector_copy(vector_mu, vector_relax_mu);
//                 vector_copy(vector_lam, vector_relax_lam);
//                 value_fun_relax = value_fun;
//                 jacobi_value_fun_relax = jacobi_value_fun;
//             }
//             flag_relax = TRUE;

//             alpha = 1.0;
//             vector_a_add_k_b(vector_x, alpha, vector_dk, vector_next_x);
//             obj_fun_next = ndscla_function_call(Objective_Fun, vector_next_x);

//             ndVectorfunction_call(con->c, vector_x, vector_temp_c);
//             for (int i = 0; i < knum_eq_constraints; i++)
//                 vector_next_h_eq->entry[i] = vector_temp_c->entry[i];
//             for (int i = 0; i < knum_ineq_constraints; i++)
//                 vector_next_g_ineq->entry[i] = vector_temp_c->entry[i + knum_eq_constraints];
//             for (int i = 0; i < knum_ineq_constraints; i++)
//                 vector_neg_g_ineq->entry[i] = max_(-1 * vector_next_g_ineq->entry[i], 0);

//             norm_next_eq_constraint = vector_1norm(vector_next_h_eq);
//             norm_next_ineq_constraint = vector_1norm(vector_next_g_ineq);

//             value_fun_next = Value_Fun(sigma, obj_fun_next, norm_next_eq_constraint, norm_next_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);
//             flag_search = TRUE;
//         }

//         // 判断下次搜索方式:
//         if (value_fun_next - value_fun_min < eta * jacobi_value_fun_min)
//         {
//             flag_search = FALSE;
//             flag_relax = FALSE;
//             relax_m = 0;
//         }

//         if (flag_relax)
//         {
//             relax_m++;
//             if (relax_m == max_relax_m && value_fun_relax < value_fun)
//             {
//                 flag_relax = FALSE;
//                 flag_search = TRUE;
//                 obj_fun = obj_fun_relax;
//                 vector_copy(vector_relax_x, vector_x);
//                 vector_copy(vector_relax_Df, vector_Df);
//                 vector_copy(vector_relax_h_eq, vector_h_eq);
//                 vector_copy(vector_relax_g_ineq, vector_g_ineq);
//                 matrix_copy(matrix_relax_Dh_eq, matrix_Dh_eq);
//                 matrix_copy(matrix_relax_Dg_ineq, matrix_Dg_ineq);
//                 matrix_copy(matrix_relax_Bk, matrix_Bk);
//                 relax_m = 0;

//                 vector_copy(vector_relax_mu, vector_mu);
//                 vector_copy(vector_relax_lam, vector_lam);

//                 value_fun = value_fun_relax;
//                 jacobi_value_fun = jacobi_value_fun_relax;
//                 continue;
//             }
//         }

//         ndscla_forward_grad(Objective_Fun, 0, vector_next_x, vector_next_Df); // TODO: 这里的参数也可能需要设置一下
//         if (value_fun_next < value_fun_min)
//         {
//             num_iteration_min = num_iterations + 1;
//             vector_copy(vector_next_x, vector_min_x);
//             value_fun_min = value_fun_next;
//             jacobi_value_fun_min = Jacobi_Value_Fun(sigma, vector_dk, vector_next_Df, norm_next_eq_constraint, norm_next_ineq_constraint, knum_eq_constraints, knum_ineq_constraints);
//             obj_fun_min = obj_fun_next;
//         }

//         // 更新自变量和拉格朗日乘子:
//         Vector *matrix_rho_dk = vector_alloc(kdimension_x);
//         vector_a_add_k_b(vector_next_x, -1, vector_x, matrix_rho_dk);

//         Flag_break = 0;
//         for (int i = 0; i < Lu; i++)
//         {
//             if (abs(matrix_rho_dk->entry[i]) >= epsilon3)
//             {
//                 Flag_break = 1;
//             }
//         }
//         if (Flag_break == 0)
//         {
//             break;
//         }

//         ndVectorfunction_jacobian(con->c, vector_x, h, matrix_temp_jacobi); // TIP 这里的差分步长需要手动设定
//         for (int i = 0; i < knum_eq_constraints; i++)
//             for (int j = 0; j < kdimension_x; j++)
//                 matrix_Dh_eq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i][j];
//         for (int i = 0; i < knum_ineq_constraints; i++)
//             for (int j = 0; j < kdimension_x; j++)
//                 matrix_Dg_ineq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i + knum_eq_constraints][j];

//         ndVectorfunction_jacobian(con->c, vector_next_x, h, matrix_temp_jacobi); // TIP 这里的差分步长需要手动设定
//         for (int i = 0; i < knum_eq_constraints; i++)
//             for (int j = 0; j < kdimension_x; j++)
//                 matrix_next_Dh_eq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i][j];
//         for (int i = 0; i < knum_ineq_constraints; i++)
//             for (int j = 0; j < kdimension_x; j++)
//                 matrix_next_Dg_ineq->matrix_entry[i][j] = matrix_temp_jacobi->matrix_entry[i + knum_eq_constraints][j];

//     }
// }