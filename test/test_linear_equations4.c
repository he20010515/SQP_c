#include <math.h>
#define N 2
#define M 2
//获取主元所在行
int getMaxinumLaber(double a[N][M], int k)
{
    int i;
    int laber = k;
    double maxinum = 0;

    //列固定，行变，找最大值
    for (i = k; i < N; i++)
        if (maxinum < fabs(a[i][k]))
        {
            maxinum = fabs(a[i][k]);
            laber = i;
        }
    return laber;
}

//列主元高斯消去法
void gaussian_elimination(double a[N][M], double *x)
{
    int i, j, k;
    int laber;
    double temp;
    // double **a = a_mai;
    // double *x=x_mai;

    //高斯消元法
    for (k = 0; k < N; k++)
    {
        laber = getMaxinumLaber(a, k);

        //主元交换,就是laber行和第k行交换
        if (laber != k)
        {
            for (i = 0; i < M; i++)
            {
                temp = a[k][i];
                a[k][i] = a[laber][i];
                a[laber][i] = temp;
            }
        }

        //消元
        for (i = k + 1; i < N; i++)
        {

            if (a[k][k] == 0)

                break;

            temp = a[i][k] / a[k][k];

            for (j = k; j < M; j++)
                a[i][j] = a[k][j] * temp - a[i][j];
        }
    }

    //回代求解x
    double sum;
    for (i = N - 1; i >= 0; i--)
    {
        sum = 0;
        for (j = i + 1; j < N; j++)
        {
            sum += a[i][j] * x[j];
        }
        x[i] = (a[i][M - 1] - sum) / a[i][i];
    }

    printf("\n");
    printf("变换增广系数矩阵\n");
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
            printf("%f   ", a[i][j]);
        printf("\n");
    }

    printf("\n");
    for (i = 0; i < N; i++)
        printf("%f   ", x[i]);
    printf("\n");
}

int main(void)
{
    double x_know[5] = {1, 2, 3, 4, 5}, y_know[5] = {7, 11, 17, 27, 40}, x[N] = {0};
    int n = 5, i = 0, j = 0;
    double A[N][M] = {0};

    for (i = 0; i < n; i++)
    {

        A[0][0] += 1;
        A[0][1] += x_know[i];
        A[0][2] += y_know[i];
        A[1][0] = A[0][1];
        A[1][1] += x_know[i] * x_know[i];
        A[1][2] += x_know[i] * y_know[i];
    }

    printf("\n");
    printf("原增广系数矩阵\n");
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
            printf("%f   ", A[i][j]);
        printf("\n");
    }
    gaussian_elimination(A, x);

    return 0;
}