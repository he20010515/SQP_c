#include "matrix.h"
#include "vector.h"

#include "elog.h"
#include "math.h"

#define LOG_TAG "simplex"

int _pivot_col(Matrix *T, int bland, double tol)
{
    // 给定一个线性规划单纯形表，确定要输入基础的变量列。

    // 参数
    // ----------
    // T : 二维数组
    //     表示单纯形表 T 的二维数组，对应于
    //     线性规划问题。它应该具有以下形式：

    //     [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
    //      [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
    //      .
    //      .
    //      .
    //      [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
    //      [c[0], c[1], ..., c[n_total], 0]]

    //     对于第 2 阶段的问题，或形式：

    //     [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
    //      [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
    //      .
    //      .
    //      .
    //      [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
    //      [c[0], c[1], ..., c[n_total], 0],
    //      [c'[0], c'[1], ..., c'[n_total], 0]]

    //      对于第 1 阶段问题（在最大化实际目标之前寻求基本可行解决方案的问题。``T`` 由``_solve_simplex`` 修改。
    // tol : 浮动
    //     目标行中大于 -tol 的元素将不考虑进行旋转。名义上这个值是零，但数值问题导致零公差是必要的。
    // Bland：布尔
    //     如果为 True，则使用 Bland 规则选择列（选择目标行中具有负系数的第一列，无论大小如何）。

    // return
    // --------
    // 状态：布尔
    //     如果找到合适的数据透视列，则为 True，否则为 False。返回 False 表示线性规划单纯形算法已完成。
    // col：整数
    //     枢轴元素的列的索引。
    //     如果 status 为 False，col 将返回为 -1;
    int flag = -1;
    int *ma = malloc(sizeof(int) * (T->col_size - 1));
    double *mma = malloc(sizeof(double) * (T->col_size - 1));
    for (int i = 0; i < T->col_size - 1; i++)
        if (T->matrix_entry[T->row_size - 1][i] < -tol)
        {
            ma[i] = true;
            mma[i] = T->matrix_entry[T->row_size - 1][i];
        }
        else
        {
            ma[i] = false;
            mma[i] = NAN;
        }
    int sum = 0;
    for (int i = 0; i < T->col_size - 1; i++)
        sum += ma[i];
    if (sum == 0)
    {
        free(ma);
        free(mma);
        return -1;
    }
    if (bland) // blad 法则,返回第一个小于-tol 的数
    {
        for (int i = 0; i < T->col_size; i++)
            if (!ma[i])
            {
                flag = i;
                break;
            }
    }
    else
    {
        Vector V;
        V.entry = mma;
        V.size = T->col_size - 1;
        flag = vector_argmin(&V);
    }
    free(ma);
    free(mma);
    return flag;
}

int _pivot_row(Matrix *T, int *basis, int basis_size, int pivcol, int phase, double tol, int bland)
{
    // 给定一个线性规划单纯形表，确定枢轴操作的行。

    // 参数
    // ----------
    // T : 二维数组
    //     一个表示单纯形表 T 的二维数组，对应于线性规划问题。它应该具有以下形式：

    //     [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
    //      [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
    //      .
    //      .
    //      .
    //      [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
    //      [c[0], c[1], ..., c[n_total], 0]]

    //     对于第 2 阶段的问题，或形式：

    //     [[A[0, 0], A[0, 1], ..., A[0, n_total], b[0]],
    //      [A[1, 0], A[1, 1], ..., A[1, n_total], b[1]],
    //      .
    //      .
    //      .
    //      [A[m, 0], A[m, 1], ..., A[m, n_total], b[m]],
    //      [c[0], c[1], ..., c[n_total], 0],
    //      [c'[0], c'[1], ..., c'[n_total], 0]]

    //      对于第一阶段问题（在最大化实际目标之前寻求基本可行解决方案的问题。“T”由“_solve_simplex”修改。
    // 基础：数组
    //     当前基本变量的列表。
    // pivcol : 整数
    //     数据透视列的索引。
    // 阶段：int
    //     单纯形算法的阶段（1 或 2）。
    // tol : 浮动
    //     枢轴列中小于 tol 的元素将不被考虑进行枢轴。名义上这个值是零，但数值问题导致零公差是必要的。
    // 平淡无奇：布尔
    //     如果为 True，则使用 Bland 规则选择行（如果可以使用多行，则选择具有最低变量索引的行）。

    // 退货
    // --------
    // 状态：布尔
    //     如果找到合适的数据透视行，则为 True，否则为 False。返回 False 表示线性规划问题是无界的。
    // 行：int
    //     枢轴元素所在行的索引。如果 status 为 False，则 row 将返回为 nan。
    // """
    int k;
    if (phase == 1)
        k = 2;
    else
        k = 1;
    int flag = -1;
    int *ma = malloc(sizeof(int) * (T->row_size - k));
    double *mma = malloc(sizeof(double) * (T->row_size - k));
    for (int i = 0; i < T->row_size - k; i++)
        if (T->matrix_entry[i][pivcol] > tol)
        {
            ma[i] = true;
            mma[i] = T->matrix_entry[i][pivcol];
        }
        else
        {
            ma[i] = false;
            mma[i] = NAN;
        }
    int sum = 0;
    for (int i = 0; i < T->row_size - k; i++)
        sum += ma[i];
    if (sum == 0)
    {
        free(ma);
        free(mma);
        return -1;
    }

    double *mmb = malloc(sizeof(double) * (T->row_size - k));
    for (int i = 0; i < T->row_size - k; i++)
    {
        if (ma[i])
            mmb[i] = T->matrix_entry[i][T->col_size - 1];
        else
            mmb[i] = NAN;
    }

    double *q = (double *)malloc(sizeof(double) * (T->row_size - k));
    for (int i = 0; i < T->row_size - k; i++)
    {
        if (ma[i])
            q[i] = mmb[i] / mma[i];
        else
            q[i] = NAN;
    }

    Vector b;
    b.size = T->row_size - k;
    b.entry = q;
    double qmin = vector_min(&b);
    int *min_rows = malloc(sizeof(int) * T->row_size - 1);
    for (int i = 0; i < T->row_size - 1; i++)
        min_rows[i] = 0;

    int len = 0;
    for (int i = 0; i < b.size; i++)
    {
        if (b.entry[i] == qmin)
        {
            min_rows[len] = i;
            len++;
        }
    }

    if (bland)
    {
        // return True, min_rows[np.argmin(np.take(basis, min_rows))]
        for (int i = 0; i < len; i++)
        {
            // if i in basis break;
            int i_in_basis = false;
            for (int j = 0; j < basis_size; j++)
            {
                if (basis[j] == i)
                {
                    i_in_basis = true;
                    break;
                }
            }
            if (i_in_basis)
            {
                flag = i;
                break;
            }
        }
    }
    else
        flag = min_rows[0];
    free(min_rows);
    free(ma);
    free(mma);
    free(mmb);
    free(q);
    return flag;
}

void _apply_pivot(Matrix *T, int *basis, int basis_size, int pivrow, int pivcol, double tol)
{
    basis[pivrow] = pivcol;
    double pivval = T->matrix_entry[pivrow][pivcol];
    for (int j = 0; j < T->col_size; j++)
        T->matrix_entry[pivrow][j] /= pivval;
    for (int irow = 0; irow < T->row_size; irow++)
    {
        if (irow != pivrow)
        {
            Vector *temp = vector_alloc(T->col_size);
            for (int i = 0; i < T->col_size; i++)
                temp->entry[i] = T->matrix_entry[pivrow][i] * T->matrix_entry[irow][pivcol];
            for (int i = 0; i < T->col_size; i++)
                T->matrix_entry[irow][i] -= temp->entry[i];
            vector_free(temp);
        }
    }
    if (fabs(pivval - tol) <= 1e-4 * tol)
    {
        log_w("The pivot operation produces a pivot value of %lf ,which is only slightly greater than the specified tolerance %lf. This may lead to issues regarding the numerical stability of the simplex method. Removing redundant constraints, changing the pivot strategy via Bland's rule or increasing the tolerance may help reduce the issue.", pivval, tol);
    }
}

int _solve_simplex(Matrix *T, int n, int *basis, int basis_size, int maxiter, double tol, int phase, int bland, int nit0, int *out_nit)
{
    int nit = nit0;
    int status = 0;
    char *message = "";
    int complete = false;
    int m = 0;
    if (phase == 1)
        m = T->col_size - 2;
    else
        m = T->col_size - 1;

    if (phase == 2)
    {
        // # Check if any artificial variables are still in the basis.
        // # If yes, check if any coefficients from this row and a column
        // # corresponding to one of the non-artificial variable is non-zero.
        // # If found, pivot at this term. If not, start phase 2.
        // # Do this for all artificial variables in the basis.
        // # Ref: "An Introduction to Linear Programming and Game Theory"
        // # by Paul R. Thie, Gerard E. Keough, 3rd Ed,
        // # Chapter 3.7 Redundant Systems (pag 102)
        int *pivrows = malloc(sizeof(int) * basis_size);
        int pivrows_len = 0;
        for (int i = 0; i < basis_size; i++)
        {
            if (basis[i] > T->col_size - 2)
            {
                pivrows[pivrows_len] = i;
                pivrows_len++;
            }
        }
        for (int i = 0; i < pivrows_len; i++) // for pivrow in temprows
        {
            int *non_zero_row = malloc(sizeof(int) * (T->col_size - 1));
            int non_zero_row_size = 0;
            for (int j = 0; j < T->col_size - 1; j++)
            {
                if (fabs(T->matrix_entry[pivrows[i]][j]) > tol)
                {
                    non_zero_row[non_zero_row_size] = j;
                    non_zero_row_size++;
                }
            }
            if (non_zero_row_size > 0)
            {
                int pivcol = non_zero_row[0];
                _apply_pivot(T, basis, basis_size, pivrows[i], pivcol, tol);
                nit++;
            }

            free(non_zero_row);
        }

        free(pivrows);
    }
    Vector *solution = NULL;
    if (m == 0)
    {
        solution = vector_alloc(T->col_size - 1);
    }
    else
    {
        int maxbasis_m = -1;
        for (int i = 0; i < basis_size; i++)
        {
            if (basis[i] >= maxbasis_m)
                maxbasis_m = basis[i];
        }
        int a = T->col_size - 1;
        int b = maxbasis_m + 1;
        solution = vector_alloc(MAX(a, b));
    }
    int pivcol = 0;
    int pivrow = 0;

    while (!complete)
    {
        pivcol = _pivot_col(T, tol, bland);
        if (pivcol == -1)
        {
            pivcol = -1;
            pivrow = -1;
            status = 0;
            complete = true;
        }
        else
        {
            pivrow = _pivot_row(T, basis, basis_size, pivcol, phase, tol, bland);
            if (pivrow == -1)
            {
                status = 3;
                complete = true;
            }
        }

        if (!complete)
        {
            if (nit >= maxiter)
            {
                // iteration limit exceeded
                status = 1;
                complete = true;
            }
            else
            {
                _apply_pivot(T, basis, basis_size, pivrow, pivcol, tol);
                nit++;
            }
        }
    }
    *out_nit = nit;
    vector_free(solution);
    return status;
}

int _linprog_simplex(const Vector *c, const Matrix *A, const Vector *b, int maxiter, double tol, int bland, Vector *x)
{
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

    // Returns
    // -------
    // x : 1-D array
    //     Solution vector.
    // status : int
    //     An integer representing the exit status of the optimization::

    //      0 : Optimization terminated successfully
    //      1 : Iteration limit reached
    //      2 : Problem appears to be infeasible
    //      3 : Problem appears to be unbounded
    //      4 : Serious numerical difficulties encountered

    // message : str
    //     A string descriptor of the exit status of the optimization.
    // iteration : int
    //     The number of iterations taken to solve the problem.

    // References
    // ----------
    // .. [1] Dantzig, George B., Linear programming and extensions. Rand
    //        Corporation Research Study Princeton Univ. Press, Princeton, NJ,
    //        1963
    // .. [2] Hillier, S.H. and Lieberman, G.J. (1995), "Introduction to
    //        Mathematical Programming", McGraw-Hill, Chapter 4.
    // .. [3] Bland, Robert G. New finite pivoting rules for the simplex method.
    //        Mathematics of Operations Research (2), 1977: pp. 103-107.

    // Notes
    // -----
    // The expected problem formulation differs between the top level ``linprog``
    // module and the method specific solvers. The method specific solvers expect a
    // problem in standard form:

    // Minimize::

    //     c @ x

    // Subject to::

    //     A @ x == b
    //         x >= 0

    // Whereas the top level ``linprog`` module expects a problem of form:

    const char *messages[] = {
        "Optimization terminated successfully.",
        "Iteration limit reached.",
        "Optimization failed. Unable to find a feasible starting point",
        "Optimization failed. The problem appears to be unbounded."
        "Optimization failed. Singular matrix encountered."};
    int n = A->row_size;
    int m = A->col_size;
    int *av = malloc(sizeof(int) * n);
    for (int i = 0; i < n; i++)
        av[i] = i + m;
    int *basis = malloc(sizeof(int) * n);

    for (int i = 0; i < n; i++)
        basis[i] = i + m;
    // fill T
    Matrix *T = matrix_alloc(n + 2, n + m + 1);
    for (int i = 0; i < T->row_size; i++)
    {
        for (int j = 0; j < T->col_size; j++)
        {
            if (i < n AND j < m) // A
                T->matrix_entry[i][j] = A->matrix_entry[i][j];
            if (i < n AND m <= j AND j < 2 * m) // E
                if (i == j - m)
                    T->matrix_entry[i][j] = 1.0;
                else
                    T->matrix_entry[i][j] = 0.0;
            if (j == T->col_size - 1 AND i < T->row_size - 1) // b
                T->matrix_entry[i][j] = b->entry[i];
            if (i == n)
                if (j < c->size)
                    T->matrix_entry[i][j] = c->entry[j];
                else
                    T->matrix_entry[i][j] = 0.0;
            if (i == n + 1)
            {
                if (j >= m AND j != T->col_size - 1)
                {
                    T->matrix_entry[i][j] = 0.0;
                }
                else
                {
                    double sum = 0.0;
                    for (int k = 0; k < n; k++)
                        sum = sum + T->matrix_entry[k][j];
                    T->matrix_entry[i][j] = -sum;
                }
            }
        }
    }
    int nit1;
    int status = _solve_simplex(T, n, basis, n, maxiter, tol, 1, bland, 0, &nit1);
    int nit2 = nit1;
    if (fabs(T->matrix_entry[T->row_size - 1][T->col_size - 1]) < tol)
    {
        // # Remove the pseudo-objective row from the tableau
        // T = T[:-1, :] //去掉最后一行
        // # Remove the artificial variable columns from the tableau
        // T = np.delete(T, av, 1) // 去掉av[0]到av[n-1]列,共n列
        Matrix *T2 = matrix_alloc(T->row_size - 1, T->col_size - n);
        for (int i = 0; i < T2->row_size; i++)
            for (int j = 0; j < T2->col_size; j++)
                if (j == T2->col_size - 1)
                    T2->matrix_entry[i][j] = T->matrix_entry[i][T->col_size - 1];

                else
                    T2->matrix_entry[i][j] = T->matrix_entry[i][j];

        matrix_free(T);
        T = T2;
    }
    else
    {
        status = 2;
        log_e("Phase 1 of the simplex method failed to find a fesible solution");
        log_e(messages[2]);
    }

    if (status == 0)
        status = _solve_simplex(T, n, basis, n, maxiter, tol, 2, bland, nit1, &nit2);
    Vector *solution = vector_alloc(n + m);
    vector_fill_const(solution, 0);
    for (int i = 0; i < n; i++)
        solution->entry[basis[i]] = T->matrix_entry[i][T->col_size - 1];
    for (int j = 0; j < m; j++)
        x->entry[j] = solution->entry[j];
    log_i("simplex optimize complete");
    log_i(messages[status]);

    free(av);
    free(basis);
    vector_free(solution);
    matrix_free(T);
    return status;
}
