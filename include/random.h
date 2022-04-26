/*
 * @Author: HeYuwei
 * @Date: 2022-03-13 09:13:40
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 18:15:34
 * @FilePath: \SQP_c\include\random.h
 * @Description:一些用于生成随机数的函数
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
/**
 * @description: 生成随机整数,左闭右开 [start,end)
 * @param {int} start
 * @param {int} end
 * @return {*}
 */
int rand_int(int start, int end);

/**
 * @description: 生成随机浮点数
 * @param {double} start
 * @param {double} end
 * @return {*}
 */
double rand_double(double start, double end);
