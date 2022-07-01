/*
 * @Author: HeYuwei
 * @Date: 2022-07-01 14:45:27
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-01 14:50:02
 * @FilePath: \SQP_c\include\sparse_matrix.h
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#pragma once

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

void s_free(storage a);
int s_get_colloc(storage a, int i);
storage s_csc2csr(storage a_csc);
void s_printf_stand(storage a);
void s_printf(storage a);