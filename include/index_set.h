/*
 * @Author: HeYuwei
 * @Date: 2022-03-13 12:36:15
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-25 11:00:12
 * @FilePath: \SQP_c\include\index_set.h
 * @Description: 指标集头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#pragma once
struct index_set
{
    int index_range;
    int *elem;
};
typedef struct index_set Index_set;

/**
 * @description:申请指标集存储空间
 * @param {int} size_指标集所能容纳的指标范围   0,1,2,...,size-1;
 * @return {*} 指向已申请指标集的指针
 */
Index_set *index_set_alloc(int size);

/**
 * @description:释放指针所指向指标集的存储空间
 * @param {Index_set} *set
 * @return {*}
 */
void index_set_free(Index_set *set);

/**
 * @description: 向指标集中添加指标
 * @param {Index_set} *set 指标集
 * @param {int} i 指标 0 <=i<= index_set.size-1
 * @return {*}
 */
void index_set_append(Index_set *set, const int i);

/**
 * @description: 在指标集中移除指标
 * @param {Index_set} *set 指标集
 * @param {int} i 指标 0 <=i<= index_set.size-1
 * @return {*}
 */
void index_set_remove(Index_set *set, const int i);

/**
 * @description: 判断指标i是否在指标集中
 * @param {Index_set} *set 指标集
 * @param {int} i 指标 0 <=i<= index_set.size-1
 * @return {*} true/false
 */
int index_set_is_in(const Index_set *set, int i);

/**
 * @description: 指标集求并集
 * @param {Index_set} *A
 * @param {Index_set} *B
 * @param {Index_set} *A_U_B
 * @return {*}
 */
void index_set_union(Index_set *A, Index_set *B, Index_set *A_U_B);

/**
 * @description: 指标集求交集
 * @param {Index_set} *A
 * @param {Index_set} *B
 * @param {Index_set} *A_I_B
 * @return {*}
 */
void index_set_intersection(const Index_set *A, const Index_set *B, Index_set *A_I_B);

/**
 * @description: 打印指标集
 * @param {Index_set} *set
 * @return {*}
 */
void index_set_print(Index_set *set);

/**
 * @description: 返回指标集当前元素的数量 !!:注意不是指标集的范围
 * @param {Index_set} *set
 * @return {*}
 */
int index_set_size(const Index_set *set);
