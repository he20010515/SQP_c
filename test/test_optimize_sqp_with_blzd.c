/*
 * @Author: HeYuwei
 * @Date: 2022-08-11 18:25:25
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-08-12 12:59:38
 * @FilePath: \SQP_c\test\test_optimize_sqp_with_blzd.c
 * @Description:
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */

#include "sqp.h"
#include "blzd_test.h"

#include "windows.h"
typedef double(__stdcall *Fun)(double[]);

int main(int argc, char const *argv[])
{
    sqp_init();

    HMODULE hlib = LoadLibraryA("D:\\Workspace\\HIT\\SQP\\BLZD\\x64\\Release\\BLZD.dll");
    if (hlib == NULL)
    {
        printf("load dll faild");
        return 0;
    }
    Fun fun = (Fun)GetProcAddress(hlib, "_target_function");
    if (fun == NULL)
    {
        printf("getProcAddress fail, error");
    }

    double u0[7] = {0};
    u0[0] = 3.8931935296046505;   //方位角
    u0[1] = 0.0267779154895057;   //极值攻角
    u0[2] = 0.0129239104473724;   //二子级俯仰角初值
    u0[3] = -0.2557582128351035;  //二子级俯仰角终端值
    u0[4] = -0.1745329251994329;  //二子级偏航角初值
    u0[5] = 0.0630259686398771;   //二子级偏航角终端值
    u0[6] = 391.3978486108092056; //二子级末端时间
    double temp;
    temp = fun(u0);
    printf("targetfunction output:%lf", temp);
    return 0;
}
