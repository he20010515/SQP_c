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
 * Print  whole entire matrix  to stdout
 *
 * @param matrix  Pointer to the matrix to be printed to stdout
 */

void matrix_print(const Matrix *matrix);

/**
 *  Printing some columns of a matrix
 *
 * @param matrix       The matrix whose part is to be printed
 * @param start_index   The  column number to start printing from
 */

void matrix_print_part(Matrix *matrix, int start_index);

/**
 * This fills the entire matrix with values which are  gotten form
 * stdin.
 *
 * @param matrix  Pointer to the matrix structure which is to
 * be filled.
 */

void matrix_fill(Matrix *matrix);

/**
 * Creates an identity matrix that is all the diagonal entries are
 * initialised to 1 and all other entries to zero).
 *
 * @param matrix_size    The number of rows and columns in the
 *                                        identity matrix.
 * @return                        A pointer to the new identity matrix just
 *                                        created
 */

Matrix *matrix_callalloc(int matrix_size);

/**
 * Allocates memory for a new matrix
 *
 * @param row_size       The number of rows in the matrix
 * @param col_size         The number of columns in the matrix
 * @return                      The location of the memory block that was
 *                                       allocated to hold the matrix
 */

Matrix *matrix_alloc(int row_size, int col_size);

/**
 * Copies the content of one matrix1 into matrix2
 *
 * @param matrix1       Pointer  the matrix to be copied
 * @param matrix2       Pointer to the matrix to which the other is to
 *                                   be copied into
 */

void matrix_copy(Matrix *matrix1, Matrix *matrix2);

/**
 * Multiplies two matrices
 *
 * @param matrix1     Pointer to the first matrix for the multiplication
 * @param matrix2     The second matrix for the multiplication
 * @return                  The location of a matrix which holds the result
 *                                of the multiplication
 */
Matrix *matrix_multiply(const Matrix *matrix1, const Matrix *matrix2);

/**
 * Multiplies a matrix by itself 'n' times
 *
 * @param matrix       The matrix which is to self multiplied 'n' times
 * @param index         The number of times the matrix is to be multiplied to
 *                                 itself
 * @return                  The location of the matrix which holds the results
 *                                    of the multiplications
 */

Matrix *matrix_pow(Matrix *matrix, int index);

/**
 * Free an entire matrix
 *
 * @param matrix       The matrix to free.
 */
void matrix_free(Matrix *matrix);

/**
 * Divides an entire row of a matrix by a value of the pivot position
 *
 * @param matrix       The matrix whose row is to be divided
 * @param pivot          The pivot position of the matrix to do the division
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

void row_operation(Matrix *multiplier_matrix, Matrix *matrix, int pivot, int row_index);

/**
 * Row echelon reduction a matrix
 *
 * @param matrix     The matrix on which to carry out row reduction
 * @param zero_control      Maximum amount of zeros that can be found on a
 *                                        row
 */

void matrix_row_reduce(Matrix *matrix, int zero_control);

/**
 * This function performs the  LU decomposition of a matrix
 *
 * @param upper_triangular        A pointer to the matrix on which to perform
 *                                                 LU decomposition.
 * @param lower_triangular          A pointer to the lower triangular matrix
 * @note You should allocate memory for the lower_triangular matrix with
 * @ref matrix_callalloc before passing it to this function
 */

void matrix_subtract(Matrix *result, Matrix *matrix1, Matrix *matrix2);

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