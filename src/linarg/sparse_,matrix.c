/*
 * @Author: HeYuwei
 * @Date: 2022-07-01 14:45:15
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-01 14:55:34
 * @FilePath: \SQP_c\src\linarg\sparse_,matrix.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include <stdlib.h>
#include <stdio.h>
#include "sparse_matrix.h"

//删除稀疏存储的矩阵   void s_free(storage a)
void s_free(storage a)
{
    free(a.type);
    free(a.index);
    free(a.indices);
    free(a.data);
}

//获取csc存储形式中indices数组第i个元素指示的元素所在的列数   int s_get_colloc(storage a, int i)
int s_get_colloc(storage a, int i)
{ //形参   第一个  结构体storage中的index数组   第二个  要查找的indices数组的第i个元素下标（即i）
    int j = 0;
    for (j; j < a.col + 1; j++)
    {
        if (*(a.index + j) > i)
            break;
    }
    return j - 1;
}

//将csc存储的数组转化为csr存储形式   storage s_csc2csr(storage a_csc)
storage s_csc2csr(storage a_csc)
{
    storage a_csr;
    a_csr.type = (char *)malloc(4 * sizeof(char));
    a_csr.type = "csr";
    a_csr.row = a_csc.row;
    a_csr.col = a_csc.col;
    a_csr.num = a_csc.num;
    a_csr.data = (double *)malloc(a_csr.num * sizeof(double));
    a_csr.indices = (int *)malloc(a_csr.num * sizeof(int));
    a_csr.index = (int *)malloc((a_csr.row + 1) * sizeof(int));

    int wmark = 0;
    a_csr.index[0] = 0;
    a_csr.index[1] = 0;
    for (int row_wmark = 0; row_wmark < a_csr.row; row_wmark++)
    {
        for (int search_rmark = 0; search_rmark < a_csc.num; search_rmark++)
        {
            if (*(a_csc.indices + search_rmark) == row_wmark)
            {
                *(a_csr.data + wmark) = *(a_csc.data + search_rmark);
                *(a_csr.index + row_wmark + 1) += 1;
                *(a_csr.indices + wmark) = s_get_colloc(a_csc, search_rmark);
                wmark++;
            }
        }
        if (row_wmark != a_csc.row - 1)
        {
            *(a_csr.index + row_wmark + 2) = *(a_csr.index + row_wmark + 1);
        }
    }
    return a_csr;
}

//输出csr存储结构体标准形式   void s_printf_stand(storage a)
void s_printf_stand(storage a)
{
    printf("====================================================\n");
    printf("data:\n");
    for (int i = 0; i < a.num; i++)
    {
        printf("%lf ", *(a.data + i));
    }
    printf("\n");
    printf("indices:\n");
    for (int i = 0; i < a.num; i++)
    {
        printf("%d ", *(a.indices + i));
    }
    printf("\n");
    printf("index:\n");
    for (int i = 0; i <= a.row; i++)
    {
        printf("%d ", *(a.index + i));
    }
    printf("\n====================================================\n");
}
//以csr格式输出存储结构体   void s_printf(storage a)
void s_printf(storage a)
{
    if (a.type == "csc")
    {
        storage a_csr = s_csc2csr(a);
        s_printf_stand(a_csr);
    }
    else
    {
        s_printf_stand(a);
    }
}

//存储结构体
typedef struct
{
    char *type;
    int num; //数据数
    int row; //行
    int col; //列
    double *data;
    int *indices; //非零元素的列索引
    int *index;   //直到第i行共有几个非零元（大小为行/列数加一，indptr[0]=0）

} storage;

//零初始化一维数组  void s_makezeros(double a[], int size)

void s_makezeros(double a[], int size)
{
    for (int i = 0; i < size; i++)
    {
        *(a + i) = 0;
    }
}

//获取二维数组中的非零元个数  int s_get_num(int row, int col, double a[][col])
int s_get_num(int row, int col, double a[][col])
{
    int num = 0;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < col; j++)
        {
            if (a[i][j] != 0)
                num += 1;
        }
    }
    return num;
}

//删除稀疏存储的矩阵   void s_free(storage a)
void s_free(storage a)
{
    free(a.type);
    free(a.index);
    free(a.indices);
    free(a.data);
}

//将二维数组转化为csr格式存储   storage s_tocsr(int row, int col, double a[][col])
storage s_tocsr(int row, int col, double a[][col])
{
    //初始化结构体
    storage s;
    s.row = row;
    s.col = col;
    s.num = s_get_num(s.row, s.col, a);
    s.type = (char *)malloc(4 * sizeof(char));
    s.type = "csr";
    s.data = (double *)malloc(s.num * sizeof(double));
    s.indices = (int *)malloc(s.num * sizeof(int));
    s.index = (int *)malloc((s.row + 1) * sizeof(int));

    int mark = 0;
    s.index[0] = 0;
    s.index[1] = 0;
    for (int i = 0; i < s.row; i++)
    {
        for (int j = 0; j < s.col; j++)
        {
            if (a[i][j] != 0)
            {
                *(s.data + mark) = a[i][j];
                *(s.indices + mark) = j;
                *(s.index + i + 1) += 1;
                mark++;
            }
        }
        if (i != s.row - 1)
        {
            *(s.index + i + 2) = *(s.index + i + 1);
        }
    }
    return s;
}

//将二维数组转化为csc格式存储    storage s_tocsc(int row, int col, double a[][col])
storage s_tocsc(int row, int col, double a[][col])
{
    storage s;
    s.row = row;
    s.col = col;
    s.num = s_get_num(s.row, s.col, a);
    s.type = (char *)malloc(4 * sizeof(char));
    s.type = "csc";
    s.data = (double *)malloc(s.num * sizeof(double));
    s.indices = (int *)malloc(s.num * sizeof(int));
    s.index = (int *)malloc((s.row + 1) * sizeof(int));

    int mark = 0;
    s.index[0] = 0;
    s.index[1] = 0;
    for (int i = 0; i < s.col; i++)
    {
        for (int j = 0; j < s.row; j++)
        {
            if (a[j][i] != 0)
            {
                *(s.data + mark) = a[j][i];
                *(s.indices + mark) = j;
                *(s.index + i + 1) += 1;
                mark++;
            }
        }
        if (i != s.col - 1)
        {
            *(s.index + i + 2) = *(s.index + i + 1);
        }
    }
    return s;
}

//获取csr存储形式中indices数组第i个元素指示的元素所在的行数   int s_get_rowloc(storage a, int i)
int s_get_rowloc(storage a, int i)
{ //形参   第一个  结构体storage中的index数组   第二个  要查找的indices数组的第i个元素下标（即i）
    int j = 0;
    for (j; j < a.row + 1; j++)
    {
        if (*(a.index + j) > i)
            break;
    }
    return j - 1;
}
//将csr存储的数组转化为csc存储形式   storage s_csr2csc(storage a_csr)
storage s_csr2csc(storage a_csr)
{
    storage a_csc;
    a_csc.type = (char *)malloc(4 * sizeof(char));
    a_csc.type = "csc";
    a_csc.row = a_csr.row;
    a_csc.col = a_csr.col;
    a_csc.num = a_csr.num;
    a_csc.data = (double *)malloc(a_csc.num * sizeof(double));
    a_csc.indices = (int *)malloc(a_csc.num * sizeof(int));
    a_csc.index = (int *)malloc((a_csc.col + 1) * sizeof(int));

    int wmark = 0;
    a_csc.index[0] = 0;
    a_csc.index[1] = 0;
    for (int col_wmark = 0; col_wmark < a_csc.col; col_wmark++)
    {
        for (int search_rmark = 0; search_rmark < a_csr.num; search_rmark++)
        {
            if (*(a_csr.indices + search_rmark) == col_wmark)
            {
                *(a_csc.data + wmark) = *(a_csr.data + search_rmark);
                *(a_csc.index + col_wmark + 1) += 1;
                *(a_csc.indices + wmark) = s_get_rowloc(a_csr, search_rmark);
                wmark++;
            }
        }
        if (col_wmark != a_csc.col - 1)
        {
            *(a_csc.index + col_wmark + 2) = *(a_csc.index + col_wmark + 1);
        }
    }
    return a_csc;
}

//获取csc存储形式中indices数组第i个元素指示的元素所在的列数   int s_get_colloc(storage a, int i)
int s_get_colloc(storage a, int i)
{ //形参   第一个  结构体storage中的index数组   第二个  要查找的indices数组的第i个元素下标（即i）
    int j = 0;
    for (j; j < a.col + 1; j++)
    {
        if (*(a.index + j) > i)
            break;
    }
    return j - 1;
}
//将csc存储的数组转化为csr存储形式   storage s_csc2csr(storage a_csc)
storage s_csc2csr(storage a_csc)
{
    storage a_csr;
    a_csr.type = (char *)malloc(4 * sizeof(char));
    a_csr.type = "csr";
    a_csr.row = a_csc.row;
    a_csr.col = a_csc.col;
    a_csr.num = a_csc.num;
    a_csr.data = (double *)malloc(a_csr.num * sizeof(double));
    a_csr.indices = (int *)malloc(a_csr.num * sizeof(int));
    a_csr.index = (int *)malloc((a_csr.row + 1) * sizeof(int));

    int wmark = 0;
    a_csr.index[0] = 0;
    a_csr.index[1] = 0;
    for (int row_wmark = 0; row_wmark < a_csr.row; row_wmark++)
    {
        for (int search_rmark = 0; search_rmark < a_csc.num; search_rmark++)
        {
            if (*(a_csc.indices + search_rmark) == row_wmark)
            {
                *(a_csr.data + wmark) = *(a_csc.data + search_rmark);
                *(a_csr.index + row_wmark + 1) += 1;
                *(a_csr.indices + wmark) = s_get_colloc(a_csc, search_rmark);
                wmark++;
            }
        }
        if (row_wmark != a_csc.row - 1)
        {
            *(a_csr.index + row_wmark + 2) = *(a_csr.index + row_wmark + 1);
        }
    }
    return a_csr;
}

//输出csr存储结构体标准形式   void s_printf_stand(storage a)
void s_printf_stand(storage a)
{
    printf("====================================================\n");
    printf("data:\n");
    for (int i = 0; i < a.num; i++)
    {
        printf("%lf ", *(a.data + i));
    }
    printf("\n");
    printf("indices:\n");
    for (int i = 0; i < a.num; i++)
    {
        printf("%d ", *(a.indices + i));
    }
    printf("\n");
    printf("index:\n");
    for (int i = 0; i <= a.row; i++)
    {
        printf("%d ", *(a.index + i));
    }
    printf("\n====================================================\n");
}
//以csr格式输出存储结构体   void s_printf(storage a)
void s_printf(storage a)
{
    if (a.type == "csc")
    {
        storage a_csr = s_csc2csr(a);
        s_printf_stand(a_csr);
    }
    else
    {
        s_printf_stand(a);
    }
}