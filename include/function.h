/*
 * @Author: HeYuwei
 * @Date: 2022-04-22 15:47:42
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-25 10:52:01
 * @FilePath: \SQP_c\include\function.h
 * @Description: 函数相关操作头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include "vector.h"
#include "matrix.h"
#pragma once
struct NdsclaFunction
{
    unsigned int inputSize;
    double (*function)(Vector *);
};
typedef struct NdsclaFunction NdsclaFunction;

/**
 * @description: 申请一个n维常量值函数
 * @param {double (*)(Vector *)} 指向一个输入Vector指针,输出Double类型的函数
 * @param {int} inputsize 函数的输入大小
 * @return {NdsclaFunction*} 指向一个NdsclaFunction对象的指针
 */
NdsclaFunction *ndscla_function_alloc(double (*function)(Vector *), int inputsize);

/**
 * @description:  调用Ndsclafunction
 * @param {NdsclaFunction} *function 函数对象
 * @param {Vector} *x x值
 * @return {double} 函数调用值
 */
double ndscla_function_call(const NdsclaFunction *function, Vector *x);

/**
 * @description: 用中心梯度法求解Ndsclafunction在x0点的梯度
 * @param {NdsclaFunction} *function
 * @param {double} h 步长
 * @param {Vector} *x0 x0点坐标
 * @param {Vector} *grad 求出梯度的存储空间
 * @return {*}
 */
void ndscla_central_grad(const NdsclaFunction *function, double h, const Vector *x0, Vector *grad);

/**
 * @description:  利用中心梯度法求解function的Hession矩阵
 * @param {NdsclaFunction} *function
 * @param {double} h 步长
 * @param {Vector} *x0 x0点坐标
 * @param {Matrix} *hession hession矩阵
 * @return {*}
 */
void ndscla_central_hession(const NdsclaFunction *function, double h, const Vector *x0, Matrix *hession);

struct vectorfunction
{
    int inputdim;  //输入维度
    int outputdim; //输出维度
    void (*function)(const Vector *, Vector *);
};
typedef struct vectorfunction NdVectorfunction;

/**
 * @description: 申请一个n维向量值函数
 * @param {void (*function)(const Vector *, Vector *)} 一个指向 具有两个Vector指针输入的函数的函数指针,其中左侧vector.size = inputdim,右侧vector.size = outputdim
 * @param {int} intputdim 输入维度
 * @param {int} outputdim 输出维度
 * @return {*}
 */
NdVectorfunction *ndVectorfunction_alloc(void (*function)(const Vector *, Vector *), int intputdim, int outputdim);

/**
 * @description: 释放一个n维向量值函数
 * @param {NdVectorfunction} *function
 * @return {*}
 */
void ndVectorfunction_free(NdVectorfunction *function);

/**
 * @description: 调用n维向量值函数
 * @param {NdVectorfunction} *function 函数对象
 * @param {Vector} *input input.size = function.inputdim
 * @param {Vector} *output output.size = function.outputdim
 * @return {*}
 */
void ndVectorfunction_call(const NdVectorfunction *function, const Vector *input, Vector *output);

/**
 * @description: 求雅克比矩阵
 * @param {NdVectorfunction} *function
 * @param {Vector} *x0
 * @param {double} h
 * @param {Matrix} *jacobian
 * @return {*}
 */
void ndVectorfunction_jacobian(const NdVectorfunction *function, const Vector *x0, double h, Matrix *jacobian);

/**
 * @description: 释放n维向量值函数存储空间
 * @param {NdsclaFunction} *function
 * @return {*}
 */
void ndscla_function_free(NdsclaFunction *function);