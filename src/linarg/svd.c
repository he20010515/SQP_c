#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char *
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1, dmaxarg2;
#define DMAX(a, b) (dmaxarg1 = (a), dmaxarg2 = (b), (dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))
static int iminarg1, iminarg2;
#define IMIN(a, b) (iminarg1 = (a), iminarg2 = (b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
    int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    double **m;
    /* allocate pointers to rows */
    m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
    m += NR_END;
    m -= nrl;
    /* allocate rows and set pointers to them */
    m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
    m[nrl] += NR_END;
    m[nrl] -= ncl;
    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;
    /* return pointer to array of pointers to rows */
    return m;
}

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
    double *v;
    v = (double *)malloc((size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
    return v - nl + NR_END;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
    free((FREE_ARG)(m[nrl] + ncl - NR_END));
    free((FREE_ARG)(m + nrl - NR_END));
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
    free((FREE_ARG)(v + nl - NR_END));
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absa > absb)
        return absa * sqrt(1.0 + (absb / absa) * (absb / absa));
    else
        return (absb == 0.0 ? 0.0 : absb * sqrt(1.0 + (absa / absb) * (absa / absb)));
}

/******************************************************************************/
void svdcmp(double **a, int m, int n, double w[], double **v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
    int flag, i, its, j, jj, k, l, nm;
    double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

    rv1 = dvector(1, n);
    g = scale = anorm = 0.0; /* Householder reduction to bidiagonal form */
    for (i = 1; i <= n; i++)
    {
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m)
        {
            for (k = i; k <= m; k++)
                scale += fabs(a[k][i]);
            if (scale)
            {
                for (k = i; k <= m; k++)
                {
                    a[k][i] /= scale;
                    s += a[k][i] * a[k][i];
                }
                f = a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = f - g;
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = i; k <= m; k++)
                        s += a[k][i] * a[k][j];
                    f = s / h;
                    for (k = i; k <= m; k++)
                        a[k][j] += f * a[k][i];
                }
                for (k = i; k <= m; k++)
                    a[k][i] *= scale;
            }
        }
        w[i] = scale * g;
        g = s = scale = 0.0;
        if (i <= m && i != n)
        {
            for (k = l; k <= n; k++)
                scale += fabs(a[i][k]);
            if (scale)
            {
                for (k = l; k <= n; k++)
                {
                    a[i][k] /= scale;
                    s += a[i][k] * a[i][k];
                }
                f = a[i][l];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][l] = f - g;
                for (k = l; k <= n; k++)
                    rv1[k] = a[i][k] / h;
                for (j = l; j <= m; j++)
                {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[j][k] * a[i][k];
                    for (k = l; k <= n; k++)
                        a[j][k] += s * rv1[k];
                }
                for (k = l; k <= n; k++)
                    a[i][k] *= scale;
            }
        }
        anorm = DMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
    }
    for (i = n; i >= 1; i--)
    { /* Accumulation of right-hand transformations. */
        if (i < n)
        {
            if (g)
            {
                for (j = l; j <= n; j++) /* Double division to avoid possible underflow. */
                    v[j][i] = (a[i][j] / a[i][l]) / g;
                for (j = l; j <= n; j++)
                {
                    for (s = 0.0, k = l; k <= n; k++)
                        s += a[i][k] * v[k][j];
                    for (k = l; k <= n; k++)
                        v[k][j] += s * v[k][i];
                }
            }
            for (j = l; j <= n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for (i = IMIN(m, n); i >= 1; i--)
    { /* Accumulation of left-hand transformations. */
        l = i + 1;
        g = w[i];
        for (j = l; j <= n; j++)
            a[i][j] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            for (j = l; j <= n; j++)
            {
                for (s = 0.0, k = l; k <= m; k++)
                    s += a[k][i] * a[k][j];
                f = (s / a[i][i]) * g;
                for (k = i; k <= m; k++)
                    a[k][j] += f * a[k][i];
            }
            for (j = i; j <= m; j++)
                a[j][i] *= g;
        }
        else
            for (j = i; j <= m; j++)
                a[j][i] = 0.0;
        ++a[i][i];
    }

    for (k = n; k >= 1; k--)
    { /* Diagonalization of the bidiagonal form. */
        for (its = 1; its <= 30; its++)
        {
            flag = 1;
            for (l = k; l >= 1; l--)
            {               /* Test for splitting. */
                nm = l - 1; /* Note that rv1[1] is always zero. */
                if ((double)(fabs(rv1[l]) + anorm) == anorm)
                {
                    flag = 0;
                    break;
                }
                if ((double)(fabs(w[nm]) + anorm) == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0; /* Cancellation of rv1[l], if l > 1. */
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    rv1[i] = c * rv1[i];
                    if ((double)(fabs(f) + anorm) == anorm)
                        break;
                    g = w[i];
                    h = pythag(f, g);
                    w[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 1; j <= m; j++)
                    {
                        y = a[j][nm];
                        z = a[j][i];
                        a[j][nm] = y * c + z * s;
                        a[j][i] = z * c - y * s;
                    }
                }
            }
            z = w[k];
            if (l == k)
            { /* Convergence. */
                if (z < 0.0)
                { /* Singular value is made nonnegative. */
                    w[k] = -z;
                    for (j = 1; j <= n; j++)
                        v[j][k] = -v[j][k];
                }
                break;
            }
            if (its == 30)
                printf("no convergence in 30 svdcmp iterations\n");
            x = w[l]; /* Shift from bottom 2-by-2 minor. */
            nm = k - 1;
            y = w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
            c = s = 1.0; /* Next QR transformation: */
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = w[i];
                h = s * g;
                g = c * g;
                z = pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 1; jj <= n; jj++)
                {
                    x = v[jj][j];
                    z = v[jj][i];
                    v[jj][j] = x * c + z * s;
                    v[jj][i] = z * c - x * s;
                }
                z = pythag(f, h);
                w[j] = z; /* Rotation can be arbitrary if z = 0. */
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 1; jj <= m; jj++)
                {
                    y = a[jj][j];
                    z = a[jj][i];
                    a[jj][j] = y * c + z * s;
                    a[jj][i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = x;
        }
    }
    free_dvector(rv1, 1, n);
}

int singular___(Matrix *mat, Matrix *matU, Matrix *matW, Matrix *matV_T)
{
    int M = mat->row_size;
    int N = mat->col_size;
    int MN;
    if (M > N)
    {
        MN = M;
    }
    else
    {
        MN = N;
    }

    double **a;
    double *w;
    double **u, **v;
    int i, j, k;
    double t;
    double *t1;
    double *t2;
    t1 = (double *)malloc(sizeof(double) * MN);
    t2 = (double *)malloc(sizeof(double) * MN);

    /* 矩阵均以M,N中的最大值来申请空间，避免越界 */
    a = dmatrix(1, MN, 1, MN);
    u = dmatrix(1, MN, 1, MN);
    w = dvector(1, MN);
    v = dmatrix(1, MN, 1, MN);

    for (int i = 0; i < mat->row_size; i++)
        for (int j = 0; j < mat->col_size; j++)
            a[i + 1][j + 1] = mat->matrix_entry[i][j];

    /*
     * 执行svdcmp后，结果U会存储在传入的矩阵a中，
     * 所以先复制矩阵a存入u中，然后在svdcmp中传入u，
     * 这样输出就保存在u中，矩阵a也不会变了
     */
    for (i = 1; i <= M; i++)
    {
        for (j = 1; j <= N; j++)
            u[i][j] = a[i][j];
    }

    // svdcmp(u, M, N, w, v);
    /*
     * 这边理应传入的是上面的参数
     * 但当 M>N时，结果不正确，U最右边的几列会全为0
     * 应该是Householder那边计算的问题，因为Householder外循环参数为N？
     *
     * 使用下面这个参数，经过验证可以得到正确的结果...
     *
     * 矩阵初始化时，元素默认赋值0，那就是说，0元素不影响Housesholder计算？后续的迭代也不影响？
     * 比如 2X2矩阵A[1,2; 2,3] 计算结果是2x2矩阵U
     * 那把矩阵改成3X3矩阵A[1,2,0; 2,3,0; 0,0,0]，即多出来的元素用0填充
     * 然后经householder计算和迭代计算后，结果U虽然是3X3矩阵，但我们还是只取2X2矩阵，两者结果是一样的？
     *
     * 用Matlab验证了下，还真的是一样的... 0元素不影响计算结果
     * 改成3x3矩阵后，U的第三列是有值的，但可以忽略掉，取前面的2x2矩阵，这样结果是一致的。
     */
    svdcmp(u, MN, MN, w, v);
    /* Sort the singular values in descending order */
    for (i = 1; i <= N; i++)
    {
        for (j = i + 1; j <= N; j++)
        {
            if (w[i] < w[j])
            { /* 对特异值排序 */
                t = w[i];
                w[i] = w[j];
                w[j] = t;
                /* 同时也要把矩阵U,V的列位置交换 */
                /* 矩阵U */
                for (k = 1; k <= M; k++)
                {
                    t1[k] = u[k][i];
                }
                for (k = 1; k <= M; k++)
                {
                    u[k][i] = u[k][j];
                }
                for (k = 1; k <= M; k++)
                {
                    u[k][j] = t1[k];
                }

                /* 矩阵V */
                for (k = 1; k <= N; k++)
                {
                    t2[k] = v[k][i];
                }
                for (k = 1; k <= N; k++)
                {
                    v[k][i] = v[k][j];
                }
                for (k = 1; k <= N; k++)
                {
                    v[k][j] = t2[k];
                }
            }
        }
    }

    /* U为MxM矩阵 */
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            matU->matrix_entry[i][j] = u[i + 1][j + 1];
        }
    }

    /* 奇异值有M个，存为W矩阵，后面会用到 */
    double **W;
    W = dmatrix(1, MN, 1, MN);
    for (i = 1; i <= M; i++)
    {
        for (j = 1; j <= N; j++)
        {
            if (i == j)
            {
                W[i][j] = w[i];
            }
            else
            {
                W[i][j] = 0.0;
            }
        }
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            matW->matrix_entry[i][j] = W[i + 1][j + 1];
        }
    }

    /* V为NxN矩阵 */
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            matV_T->matrix_entry[i][j] = v[j + 1][i + 1];
        }
    }

    /* 验证结果，即计算U*W*V’看是否等于a */

    /* 再算 temp*V */
    /* 先对v进行矩阵转置，存为V */

    free_dmatrix(a, 1, MN, 1, MN);
    free_dmatrix(u, 1, MN, 1, MN);
    free_dvector(w, 1, MN);
    free_dmatrix(v, 1, MN, 1, MN);
    free_dmatrix(W, 1, MN, 1, MN);
    free(t1);
    free(t2);

    return 0;
}

void matrix_singluar(const Matrix *A, Matrix *U, Matrix *W, Matrix *V_T)
{
    // 矩阵的奇异值分解
    // A m*n
    // U m*m
    // W m*n
    // V_T n*n
    int m = A->row_size;
    int n = A->col_size;
    if (!((m == U->row_size AND m == U->col_size)
              AND(m == W->row_size AND n == W->col_size)
                  AND(n == V_T->row_size AND n == V_T->col_size)))
    {
        terminate("SVD ERROR");
    }
    singular___(A, U, W, V_T);
}
