#include <stdio.h>
int main()
{
    double A[4][4];
    double L[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    double U[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    double Lni[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    double Uni[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    double Ani[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
    int size = 4;

    printf("输入四阶矩阵A:\n");
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            scanf("%lf", &A[i][j]);

    for (int i = 0; i < size; i++)
    {
        L[i][i] = 1.0;
    }
    for (int j = 0; j < size; j++)
    {
        U[0][j] = A[0][j];
    }
    for (int i = 1; i < size; i++)
    {
        L[i][0] = A[i][0] / U[0][0];
    }
    for (int k = 1; k < size; k++)
    {
        for (int j = k; j < size; j++)
        {
            double s = 0.0;
            for (int t = 0; t < k; t++)
            {
                s += L[k][t] * U[t][j];
            }
            U[k][j] = A[k][j] - s;
        }
        for (int i = k; i < size; i++)
        {
            double s = 0.0;
            for (int t = 0; t < k; t++)
            {
                s += L[i][t] * U[t][k];
            }
            L[i][k] = (A[i][k] - s) / U[k][k];
        }
    }
    for (int j = 0; j < size; j++)
    {
        for (int i = j; i < size; i++)
        {
            if (i == j)
                Lni[i][j] = 1 / L[i][j];
            else if (i < j)
                Lni[i][j] = 0;
            else
            {
                double s = 0.0;
                for (int k = j; k < i; k++)
                {
                    s += L[i][k] * Lni[k][j];
                }
                Lni[i][j] = -Lni[j][j] * s;
            }
        }
    }
    for (int j = 0; j < size; j++)
    {
        for (int i = j; i >= 0; i--)
        {
            if (i == j)
                Uni[i][j] = 1 / U[i][j];
            else if (i > j)
                Uni[i][j] = 0;
            else
            {
                double s = 0.0;
                for (int k = i + 1; k <= j; k++)
                {
                    s += U[i][k] * Uni[k][j];
                }
                Uni[i][j] = -1 / U[i][i] * s;
            }
        }
    }
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            for (int k = 0; k < size; k++)
            {
                Ani[i][j] += Uni[i][k] * Lni[k][j];
            }
        }
    }

    printf("L矩阵:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf(" %8.2f", L[i][j]);
        }
        printf("\n");
    }
    printf("U矩阵:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf(" %8.2f", U[i][j]);
        }
        printf("\n");
    }
    printf("L的逆矩阵:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf(" %8.2f", Lni[i][j]);
        }
        printf("\n");
    }
    printf("U的逆矩阵:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf(" %8.2f", Uni[i][j]);
        }
        printf("\n");
    }
    printf("A的逆矩阵:\n");
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf(" %8.2f", Ani[i][j]);
        }
        printf("\n");
    }
    return 0;
}
