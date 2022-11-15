/*
 * @Author: heyuwei he20010515@163.com
 * @Date: 2022-09-18 18:59:52
 * @LastEditors: heyuwei he20010515@163.com
 * @LastEditTime: 2022-11-15 12:32:50
 * @FilePath: /SQP_c/src/util/random.c
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
#include <math.h>
int rand_int(int start, int end)
{
    // 返回一个 start 到 end的随机数,左闭右开 [start,end)
    if (end < start)
    {
        terminate("ERROR rand_int,start must less than < end");
    }
    int range = end - start;
    int k = rand() % range; // 计算整数k
    return k + start;
}

double rand_double(double start, double end)
{
    // 返回一个 start 到 end 的随机double浮点数,左闭右开 [start,end)
    if (end < start)
    {
        terminate("ERROR rand_int,start must less than < end");
    }
    double x = (double)rand() / RAND_MAX * (end - start) + start;
    return x;
}

#define PI 3.1415926535

/**
 * @description: 产生标准正态分布 均值为0 标准差为 1
 * @return {*}
 */
double rand_gauss()

{
    static double U, V;
    static int phase = 0;
    double z;
    U = rand() / (RAND_MAX + 1.0);
    V = rand() / (RAND_MAX + 1.0);
    z = sqrt(-2.0 * log(U)) * sin(2.0 * PI * V);

    return z;
}