/*
 * @Author: HeYuwei
 * @Date: 2022-03-13 09:13:40
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-07-08 11:22:24
 * @FilePath: \SQP_c\include\util.h
 * @Description:一些常用工具函数和预定义宏
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "vector.h"
#pragma once

#define AND &&
#define OR ||
#define DOUBLE_ERROR 0.01
/**
 * @description:出现错误退出的函数,引发异常可以中断并调试
 * @param {char} *string
 * @return {*}
 */
void terminate(char *string);

/**
 * @description: 判断浮点数是否相等
 * @param {double} a
 * @param {double} b
 * @return {*}
 */
int double_equal(double a, double b);

/**
 * @description: 判断向量中元素是否都大于一个常数a
 * @param {Vector} *v
 * @param {double} a
 * @return {*}
 */
int vector_any_bigger_equal_than_const(const Vector *v, double a);

void opem_mp_test(void);

/**
 * @description: 初始化函数,用于初始化日志系统
 * @param {*}
 * @return {*}
 */
void sqp_init(void);

#define TRUE 1
#define FALSE 0