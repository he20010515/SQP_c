/*
 * @Author: HeYuwei
 * @Date: 2022-03-27 19:10:22
 * @LastEditors: Heyuwei
 * @LastEditTime: 2022-04-26 18:09:02
 * @FilePath: \SQP_c\include\matrix.h
 * @Description: 矩阵/向量联合运算头文件
 *
 * Copyright (c) 2022 by Heyuwei, All Rights Reserved.
 */
#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "vector.h"
#include "index_set.h"

#pragma once

/**
 * Definition of a  @ref _Matrix structure to hold the matrix entries.
 * A new matrix can be created using the @ref matrix_alloc function.
 */
struct _Matrix
{
  /** The number  of rows in the matrix */
  int row_size;

  /** The number of columns in the matrix */
  int col_size;

  /** The individual entries in the matrix */
  double **matrix_entry;
};

/**
 * Represents the @link _Matrix matrix structure@endlink which holds the matrix entries.
 * to initialise a new Matrix, simply create a variable of this type
 * and initialise to NULL if necessary
 *
 * @see matrix_alloc
 */

typedef struct _Matrix Matrix;

/**
 * @description: 打印matrix所指向的矩阵
 * @param {Matrix} *matrix
 * @return {*}
 */
void matrix_print(const Matrix *matrix);

/**
 * @description: 打印matrix的一部分
 * @param {Matrix} *matrix
 * @param {int} start_index
 * @return {*}
 */
void matrix_print_part(Matrix *matrix, int start_index);

/**
 * @description: 通过scanf手动填充矩阵
 * @param {Matrix} *matrix
 * @return {*}
 */
void matrix_fill(Matrix *matrix);

/**
 * @description: 申请一个单位阵
 * @param {int} matrix_size
 * @return {*}
 */
Matrix *matrix_callalloc(int matrix_size);

/**
 * @description: 申请一个矩阵并分配存储空间
 * @param {int} row_size 行
 * @param {int} col_size 列
 * @return {*}
 */
Matrix *matrix_alloc(int row_size, int col_size);

/**
 * @description: 矩阵复制,将matrix1复制到matrix2
 * @param {Matrix} *matrix1 源 矩阵指针
 * @param {Matrix} *matrix2 目标矩阵指针
 * @return {*}
 */
void matrix_copy(Matrix *matrix1, Matrix *matrix2);

/**
 * @description: 矩阵乘法,将matrix1与matrix2相乘,返回一个指向结果的指针
 * @param {Matrix} *matrix1
 * @param {Matrix} *matrix2
 * @return {*} 乘法运算结果,需要调用者手动释放
 */
Matrix *matrix_multiply(const Matrix *matrix1, const Matrix *matrix2);

/**
 * @description: 矩阵求幂
 * @param {Matrix} *matrix
 * @param {int} index
 * @return {*} 求幂运算的结果,需要调用者手动释放
 */
Matrix *matrix_pow(Matrix *matrix, int index);

/**
 * @description: 释放存储空间
 * @param {Matrix} *matrix
 * @return {*}
 */
void matrix_free(Matrix *matrix);

/**
 * @description: 矩阵的第i行除以第i行,第pivot列的值
 * @param {Matrix} *matrix
 * @param {int} pivot 列坐标
 * @return {*}
 */
void row_divide(Matrix *matrix, int pivot);

/**
 * Row operations on the matrix
 *
 * @param multiplier_matrix    A matrix to store the various multipliers used
 *                                             in row reduction
 * @param matrix          A matrix on which to carry out the row
 *                                  operations
 * @param pivot            The pivot position of the matrix to use
 * @param  row_index   The row number on which to carry out row operations
 */

/**
 * @description: 对矩阵的行进行操作
 * @param {Matrix} *multiplier_matrix
 * @param {Matrix} *matrix
 * @param {int} pivot
 * @param {int} row_index
 * @return {*}
 */
void row_operation(Matrix *multiplier_matrix, Matrix *matrix, int pivot, int row_index);

/**
 * Adds one matrix to another
 *
 * @param result        A  matrix to hold the result of the addition
 * @param matrix1     The first matrix for the addition
 * @param matrix2     The second matrix for the addition
 */

void matrix_add(Matrix *result, Matrix *matrix1, Matrix *matrix2);

/**
 * Checks if two matrices have equal rows and columns
 *
 * @param matrix1    The first matrix
 * @param matrix2    The second matrix
 * @return          Non-zero if the matrix row and columns  are equal, zero if they
 *                        are not equal.
 */
int matrix_equal_size(Matrix *matrix1, Matrix *matrix2);

/**
 * Checks  if there are too many zeros in a single line
 *
 * @param matrix The matrix which is to be checked
 * @param control_index    The maximum amount of zero's that can be contained
 *                                         in a single row
 */
void error_zeros(Matrix *matrix, int control_index);

/**
 * Function to terminate an application in case of an error
 *
 *  @param string     Message to displayed to stdout in case of an error
 */
void terminate(char *string);

/**
 * Function to compute F2-norm of A
 *
 *  @param matrix     A Matrix
 *
 *  @return double    F2norm of a matrix
 */
double matrix_F2norm(Matrix *matrix);

void matrix_inverse(Matrix *mat, Matrix *inv);

void array_2_matrix(double *array, const int rowsize, const int colsize, Matrix *mat);

void matrix_fill_const(Matrix *mat, double a);

void matrix_lu_depose(Matrix *mat, Matrix *L, Matrix *U);

void matrix_mutiply_vector(const Matrix *mat, const Vector *a, Vector *mat_a);

void vector_mutiply_matrix(const Vector *a, const Matrix *mat, Vector *mat_a);

void matrix_submatrix_by_rowindex_set(const Matrix *A, const Index_set *index_set, Matrix *subA);

void matrix_transpose(const Matrix *mat, Matrix *matT);

void vector_mutiply_vectorT(const Vector *V, const Vector *W, Matrix *VWT);

void matrix_mutiply_const(const Matrix *A, double k, Matrix *kA);

void vector_log(const Vector *v);
void matrix_log(const Matrix *mat);

int matrix_have_na(const Matrix *mat);
int vector_have_na(const Vector *v);