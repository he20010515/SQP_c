#include <math.h>
#include "string.h"
#include "matrix.h"
#include "util.h"
#include "vector.h"
#include "index_set.h"
void matrix_print(const Matrix *matrix)
{
	int i, j;
	printf("\n");
	for (i = 0; i < matrix->row_size; i++)
	{
		printf("\t\t");
		for (j = 0; j < matrix->col_size; j++)
		{
			printf("%9.2g", matrix->matrix_entry[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

void matrix_print_part(Matrix *matrix, int start_index)
{
	int j, i;
	for (i = 0; i < matrix->row_size; ++i)
	{
		for (j = start_index; j < matrix->col_size; ++j)
		{
			printf("\t\t%9.2f", matrix->matrix_entry[i][j]);
		}
		printf("\n");
	}
}

void matrix_fill(Matrix *matrix)
{
	int i, j;
	printf("Enter the contents of the matrix:\n");
	for (i = 0; i < matrix->row_size; i++)
	{
		for (j = 0; j < matrix->col_size; j++)
		{
			scanf("%lf", &matrix->matrix_entry[i][j]);
		}
	}
}

/*Function to create an identity matrix	*/
Matrix *matrix_callalloc(int matrix_size)
{
	Matrix *result = matrix_alloc(matrix_size, matrix_size);
	int i, j;

	for (i = 0; i < matrix_size; i += 1)
	{
		for (j = 0; j < matrix_size; j += 1)
		{
			if (j == i)
			{
				result->matrix_entry[i][j] = 1;
			}

			else
			{
				result->matrix_entry[i][j] = 0;
			}
		}
	}

	return result;
}

Matrix *matrix_alloc(int row_size, int col_size)
{
	int j;
	Matrix *new_matrix = malloc(sizeof(Matrix));

	// Allocating memory for the new matrix structure
	new_matrix->row_size = row_size;
	new_matrix->col_size = col_size;
	new_matrix->matrix_entry = malloc(new_matrix->row_size * sizeof(double *));
	for (j = 0; j < new_matrix->row_size; j++)
	{
		new_matrix->matrix_entry[j] = malloc(new_matrix->col_size * sizeof(double));
	}
	for (int i = 0; i < row_size; i++)
	{
		for (int j = 0; j < col_size; j++)
		{
			new_matrix->matrix_entry[i][j] == 0.0;
		}
	}

	return new_matrix;
}

/*Copies Matrix1 into matrix2 */
void matrix_copy(Matrix *matrix1, Matrix *matrix2)
{
	int i, j;
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			matrix2->matrix_entry[i][j] = matrix1->matrix_entry[i][j];
		}
	}
}

Matrix *matrix_multiply(const Matrix *matrix1, const Matrix *matrix2)
{
	int i, j, k;
	double sum;
	if (matrix1->col_size != matrix2->row_size)
	{
		terminate("ERROR: The number columns of matrix1  != number of rows in matrix2!");
	}
	Matrix *result = matrix_alloc(matrix1->row_size, matrix2->col_size);
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (k = 0; k < matrix2->col_size; k += 1)
		{
			sum = 0.0;

			for (j = 0; j < matrix1->col_size; j += 1)
			{
				sum += matrix1->matrix_entry[i][j] * matrix2->matrix_entry[j][k];
			}

			result->matrix_entry[i][k] = sum;
		}
	}
	return result;
}

Matrix *matrix_pow(Matrix *matrix, int index)
{
	if (index == 1)
	{
		Matrix *result = matrix_alloc(matrix->row_size, matrix->col_size);
		matrix_copy(matrix, result);
		return result;
	}
	else
	{

		int i, j, k, l, sum, count;

		Matrix *temp = matrix_alloc(matrix->row_size, matrix->col_size);   // Allocating space for a temporal matrix
		Matrix *result = matrix_alloc(matrix->row_size, matrix->col_size); // Allocating space for the result matrix

		matrix_copy(matrix, temp);

		count = index / 2 - 1;
		if (count < 1)
		{
			matrix_copy(matrix, result);
		}

		else
		{
			for (l = 0; l < count; l += 1)
			{
				for (i = 0; i < matrix->row_size; i += 1)
				{
					for (k = 0; k < matrix->col_size; k += 1)
					{
						sum = 0;

						for (j = 0; j < matrix->col_size; j += 1)
						{
							sum += (temp->matrix_entry[i][j] * matrix->matrix_entry[j][k]);
						}

						result->matrix_entry[i][k] = sum;
					}
				}

				/* Copying the result matrix into the temp matrix for further
				 * multiplication */
				matrix_copy(result, temp);
			}
		}

		/*	Freeing the temp matrix		*/
		matrix_free(temp);
		if (index % 2 == 0)
		{
			Matrix *result_final = matrix_multiply(result, result);
			/* Freeing the result Matrix	*/
			matrix_free(result);

			return result_final;
		}

		else
		{
			Matrix *temp = matrix_multiply(matrix, result);
			Matrix *result_final = matrix_multiply(temp, result);

			/* Freeing the temp matrix		*/
			matrix_free(temp);

			/* Freeing the result Matrix	*/
			matrix_free(result);

			return result_final;
		} // End of else statement
	}
}

void matrix_free(Matrix *matrix)
{
	int j;
	for (j = 0; j < matrix->row_size; j++)
	{
		free(matrix->matrix_entry[j]);
	}
	free(matrix->matrix_entry);
	free(matrix);
}

/*Function which divides all row entries by the value of a the diagonal */
void row_divide(Matrix *matrix, int pivot)
{
	int j;
	double divisor = matrix->matrix_entry[pivot][pivot],
		   result;

	for (j = pivot; j < matrix->col_size; j++)
	{
		result = (matrix->matrix_entry[pivot][j] / divisor);
		matrix->matrix_entry[pivot][j] = result;
	}
}

/*Function to carry out row operations*/
void row_operation(Matrix *multiplier_matrix, Matrix *matrix, int pivot, int row_index)
{
	int j;
	double multiplier = (matrix->matrix_entry[row_index][pivot] / matrix->matrix_entry[pivot][pivot]);
	// Loop which checks if matrix is provided to store the multiplier
	if (multiplier_matrix != NULL)
	{
		multiplier_matrix->matrix_entry[row_index][pivot] = multiplier;
	}

	for (j = 0; j < matrix->col_size; j++)
	{
		matrix->matrix_entry[row_index][j] -= multiplier * matrix->matrix_entry[pivot][j];
	}
}

void matrix_row_reduce(Matrix *matrix, int zero_control)
{
	int pivot, row_index;
	double multiplier;
	for (pivot = 0; pivot < matrix->row_size; pivot++)
	{

		error_zeros(matrix, zero_control); // Function checks if there are too many zeros in a single row
		if ((matrix->matrix_entry[pivot][pivot] != 1) && (matrix->matrix_entry[pivot][pivot] != 0))
		{
			row_divide(matrix, pivot);
		}

		for (row_index = pivot + 1; row_index < matrix->row_size; row_index++)
		{
			if (matrix->matrix_entry[pivot][pivot] != 0)
			{
				row_operation(NULL, matrix, pivot, row_index);
			}
		}

		for (row_index = pivot - 1; row_index >= 0; row_index--)
		{
			if (matrix->matrix_entry[pivot][pivot] != 0)
			{
				row_operation(NULL, matrix, pivot, row_index);
			}
		}
	}
}

void matrix_subtract(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
	int i, j;

	if (!(matrix_equal_size(matrix1, matrix2)) ||
		!(matrix_equal_size(matrix2, result)))
	{
		terminate("ERROR: The matrices you are trying to subtract have different sizes");
	}

	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] - matrix2->matrix_entry[i][j];
		}
	}
}

void matrix_add(Matrix *result, Matrix *matrix1, Matrix *matrix2)
{
	int i, j;
	if (!(matrix_equal_size(matrix1, matrix2)) ||
		!(matrix_equal_size(matrix2, result)))
	{
		terminate("ERROR: The matrices you are trying to add  have different sizes");
	}
	for (i = 0; i < matrix1->row_size; i += 1)
	{
		for (j = 0; j < matrix1->col_size; j += 1)
		{
			result->matrix_entry[i][j] = matrix1->matrix_entry[i][j] + matrix2->matrix_entry[i][j];
		}
	}
}

void matrix_inverse(Matrix *mat, Matrix *inv)
{
	if (mat->col_size != mat->row_size)
	{
		terminate("ERROR: The matrix to inverse must have same colsize and rowsize");
	}
	int size = mat->row_size;
	Matrix *L = matrix_alloc(size, size);
	matrix_fill_const(L, 0.0);
	Matrix *U = matrix_alloc(size, size);
	matrix_fill_const(U, 0.0);
	Matrix *Lni = matrix_alloc(size, size);
	matrix_fill_const(Lni, 0.0);
	Matrix *Uni = matrix_alloc(size, size);
	matrix_fill_const(Uni, 0.0);

	double s;
	for (int i = 0; i < size; i++)
	{
		L->matrix_entry[i][i] = 1.0;
	}
	for (int j = 0; j < size; j++)
	{
		U->matrix_entry[0][j] = mat->matrix_entry[0][j];
	}
	for (int i = 1; i < size; i++)
	{
		L->matrix_entry[i][0] = mat->matrix_entry[i][0] / U->matrix_entry[0][0];
	}
	for (int k = 1; k < size; k++)
	{
		for (int j = k; j < size; j++)
		{
			s = 0.0;
			for (int t = 0; t < k; t++)
			{
				s += L->matrix_entry[k][t] * U->matrix_entry[t][j];
			}
			U->matrix_entry[k][j] = mat->matrix_entry[k][j] - s;
		}
		for (int i = k; i < size; i++)
		{
			s = 0.0;
			for (int t = 0; t < k; t++)
			{
				s += L->matrix_entry[i][t] * U->matrix_entry[t][k];
			}
			L->matrix_entry[i][k] = (mat->matrix_entry[i][k] - s) / U->matrix_entry[k][k];
		}
	}
	for (int j = 0; j < size; j++)
	{
		for (int i = j; i < size; i++)
		{
			if (i == j)
				Lni->matrix_entry[i][j] = 1. / L->matrix_entry[i][j];
			else if (i < j)
				Lni->matrix_entry[i][j] = 0.;
			else
			{
				s = 0.0;
				for (int k = j; k < i; k++)
				{
					s += L->matrix_entry[i][k] * Lni->matrix_entry[k][j];
				}
				Lni->matrix_entry[i][j] = -Lni->matrix_entry[j][j] * s;
			}
		}
	}
	for (int j = 0; j < size; j++)
	{
		for (int i = j; i >= 0; i--)
		{
			if (i == j)
				Uni->matrix_entry[i][j] = 1. / U->matrix_entry[i][j];
			else if (i > j)
				Uni->matrix_entry[i][j] = 0.;
			else
			{
				s = 0.0;
				for (int k = i + 1; k <= j; k++)
				{
					s += U->matrix_entry[i][k] * Uni->matrix_entry[k][j];
				}
				Uni->matrix_entry[i][j] = -1. / U->matrix_entry[i][i] * s;
			}
		}
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				inv->matrix_entry[i][j] += Uni->matrix_entry[i][k] * Lni->matrix_entry[k][j];
			}
		}
	}
	matrix_free(L);
	matrix_free(U);
	matrix_free(Lni);
	matrix_free(Uni);
	return;
}

int matrix_equal_size(Matrix *matrix1, Matrix *matrix2)
{

	return (matrix1->row_size == matrix2->row_size &&
			matrix1->col_size == matrix2->col_size);
}

/*
  This function checks if there is a line containing too many zero's and it exits
  if such a line is found
*/
void error_zeros(Matrix *matrix, int control_index)
{
	int i, j, count;
	for (i = 0; i < matrix->row_size; i++)
	{
		count = 0;
		for (j = 0; j < matrix->col_size; j++)
		{
			if (matrix->matrix_entry[i][j] == 0)
			{
				count++;
			}
			else
			{
				return;
			}
			if (count == control_index)
			{
				fprintf(stdout, "\nProcess fail because row %d contains %d  zeros\n", i + 1, control_index);
				matrix_print(matrix);
				exit(1);
			}
		}
	}
}
double matrix_F2norm(Matrix *matrix)
{
	double norm = 0.;
	for (int i = 0; i < matrix->row_size; i++)
	{
		for (int j = 0; j < matrix->col_size; j++)
		{
			norm += pow(matrix->matrix_entry[i][j], 2.);
		}
	}
	return sqrt(norm);
}

void array_2_matrix(double *array, const int rowsize, const int colsize, Matrix *mat)
{
	if ((rowsize != mat->row_size) || (colsize != mat->col_size))
	{
		terminate("ERROR array_2_matrix");
	}

	for (size_t i = 0; i < rowsize; i++)
	{
		for (size_t j = 0; j < colsize; j++)
		{
			mat->matrix_entry[i][j] = array[colsize * i + j];
		}
	}
	return;
}

void matrix_fill_const(Matrix *mat, double a)
{
	for (size_t i = 0; i < mat->row_size; i++)
	{
		for (size_t j = 0; j < mat->col_size; j++)
		{
			mat->matrix_entry[i][j] = a;
		}
	}
}

void matrix_lu_depose(Matrix *mat, Matrix *L, Matrix *U)
{
	if (!(matrix_equal_size(mat, L) AND matrix_equal_size(L, U)))
	{
		terminate("ERROR LUdepose must have same size");
	}

	int size = mat->row_size;
	matrix_fill_const(L, 0.0);
	matrix_fill_const(U, 0.0);

	double s;
	for (int i = 0; i < size; i++)
	{
		L->matrix_entry[i][i] = 1.0;
	}
	for (int j = 0; j < size; j++)
	{
		U->matrix_entry[0][j] = mat->matrix_entry[0][j];
	}
	for (int i = 1; i < size; i++)
	{
		L->matrix_entry[i][0] = mat->matrix_entry[i][0] / U->matrix_entry[0][0];
	}
	for (int k = 1; k < size; k++)
	{
		for (int j = k; j < size; j++)
		{
			s = 0.0;
			for (int t = 0; t < k; t++)
			{
				s += L->matrix_entry[k][t] * U->matrix_entry[t][j];
			}
			U->matrix_entry[k][j] = mat->matrix_entry[k][j] - s;
		}
		for (int i = k; i < size; i++)
		{
			s = 0.0;
			for (int t = 0; t < k; t++)
			{
				s += L->matrix_entry[i][t] * U->matrix_entry[t][k];
			}
			L->matrix_entry[i][k] = (mat->matrix_entry[i][k] - s) / U->matrix_entry[k][k];
		}
	}

	return;
}

void matrix_mutiply_vector(const Matrix *mat, const Vector *a, Vector *mat_a)
{
	if (mat->col_size != a->size)
	{
		terminate("ERROR matrix_mutiply_vector: the matrix.col must equal to vector.size");
	}
	if (mat_a->size != mat->row_size)
	{
		terminate("ERROR matrix_mutiply_vector: a.size != mat_a.size ");
	}

	double s = 0.0;
	for (size_t i = 0; i < mat_a->size; i++)
	{
		s = 0.0;
		for (size_t j = 0; j < mat->col_size; j++)
		{
			s = s + mat->matrix_entry[i][j] * a->entry[j];
		}
		mat_a->entry[i] = s;
	}
}

void vector_mutiply_matrix(Vector *a, Matrix *mat, Vector *mat_a)
{
	if (mat->col_size != a->size)
	{
		terminate("ERROR matrix_mutiply_vector: the matrix.col_size must equal to vector.size");
	}
	if (a->size != mat_a->size)
	{
		terminate("ERROR matrix_mutiply_vector: a.size != mat_a.size ");
	}

	double s = 0.0;
	for (size_t j = 0; j < a->size; j++)
	{
		s = 0.0;
		for (size_t i = 0; i < mat->col_size; i++)
		{
			s = s + mat->matrix_entry[i][j] * a->entry[i];
		}
		mat_a->entry[j] = s;
	}
}

void matrix_submatrix_by_rowindex_set(const Matrix *A, const Index_set *index_set, Matrix *subA)
{
	if (!(A->col_size == subA->col_size AND index_set->index_range == A->row_size))
	{
		terminate("ERROR matrix_submatrix_by_rouindex_set size not fit");
	}
	int sub_i = 0;
	for (int i = 0; i < A->row_size; i++)
	{
		if (index_set_is_in(index_set, i))
		{
			for (int j = 0; j < A->col_size; j++)
			{
				subA->matrix_entry[sub_i][j] = A->matrix_entry[i][j];
			}
			sub_i++;
		}
	}
}

void matrix_transpose(const Matrix *mat, Matrix *matT)
{
	if (!(mat->col_size == matT->row_size AND mat->row_size == matT->col_size))
	{
		terminate("ERROR matrix_transpose matrix size not fit");
	}
	for (int i = 0; i < mat->row_size; i++)
	{
		for (int j = 0; j < mat->col_size; j++)
		{
			matT->matrix_entry[j][i] = mat->matrix_entry[i][j];
		}
	}
}