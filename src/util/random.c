#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "util.h"
int rand_int(int start, int end)
{
    //返回一个 start 到 end的随机数,左闭右开 [start,end)
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
    //返回一个 start 到 end 的随机double浮点数,左闭右开 [start,end)
    if (end < start)
    {
        terminate("ERROR rand_int,start must less than < end");
    }
    double x = (double)rand() / RAND_MAX * (end - start) + start;
    return x;
}