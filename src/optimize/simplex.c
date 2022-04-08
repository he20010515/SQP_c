#include "matrix.h"
#include "vector.h"
#include "elog.h"
#include "lp.h"
#include "math.h"
#include "simplex.h"

#define LOG_TAG "simplex"

/*
Solve Problem :
    min         z = c^T x
    s.t.   A_eq x = b_eq
           A_ub x = b_ub
           x>=0,b_eq>=0,c_eq>=0;
Method:
Simplex *problem = simplex_alloc(c,A_ub,b_ub,A_eq,B_eq);
int flag = simplex_main(problem,xstar);
    //* flag 1 Optimal solution found
    //* flag 0 infity solution
    //*
simplex_free(problem)
*/

Simplex *simplex_alloc(const Vector *c, const Matrix *A_ub, const Vector *b_ub, const Matrix *A_eq, const Vector *b_eq)
{
    Simplex *self = (Simplex *)malloc(sizeof(Simplex));
    self->c = vector_alloc(c->size);
    vector_copy(c, self->c);
    if (A_ub == NULL OR b_ub == NULL)
    {
        self->A_ub = NULL;
        self->b_ub = NULL;
    }
    else
    {
        self->A_ub = matrix_alloc(A_ub->row_size, A_ub->col_size);
        matrix_copy((Matrix *)A_ub, self->A_ub);
        self->b_ub = vector_alloc(b_ub->size);
        vector_copy(b_ub, self->b_ub);
        assert((A_ub->row_size == b_ub->size) AND(A_ub->col_size == c->size))
    }
    if (A_eq == NULL OR b_eq == NULL)
    {
        self->A_eq = NULL;
        self->b_eq = NULL;
    }
    else
    {
        self->A_eq = matrix_alloc(A_eq->row_size, A_eq->col_size);
        matrix_copy((Matrix *)A_eq, self->A_eq);
        self->b_eq = vector_alloc(b_eq->size);
        vector_copy(b_eq, self->b_eq);
        assert((A_eq->row_size == b_eq->size) AND(A_eq->col_size == c->size))
    }

    if (self->A_eq == NULL AND self->A_ub == NULL)
    {
        log_e("must have linearconstraints");
        exit(-1);
    }

    self->N_x = 0;
    self->meq = 0;
    self->mub = 0;
    self->m = 0;
    self->T = NULL;
    self->sign = NULL;
    self->F = 1;
}

void simplex_inital_value(Simplex *self)
{
    self->N_x = self->c->size;
    // 按照约束类型求解
    if (self->A_eq != NULL) //判断等式约束有无
        self->meq = self->A_eq->row_size;
    if (self->A_ub != NULL) //判断不等式约束有无
        self->mub = self->A_ub->row_size;
    self->m = self->mub + self->meq; //所有约束的行数
    // 建立第一阶段单纯形表
    self->T = matrix_alloc(self->m + 1, self->N_x + self->m + 1);
    matrix_fill_const(self->T, 0.0);
    // 判断两种约束的有无
    // T矩阵上半部分为等式约束，下半部分为不等式约束
    // b = self.T[:-1,-1]
    if (self->meq > 0)
    {
        // self.T[:self.meq, :self.N_x] = self.A_eq
        for (int i = 0; i < self->meq; i++)
            for (int j = 0; j < self->N_x; j++)
                self->T->matrix_entry[i][j] = self->A_eq->matrix_entry[i][j];
        // b[:self.meq] = self.b_eq
        for (int i = 0; i < self->meq; i++)
            self->T->matrix_entry[i][self->T->col_size - 1] = self->b_eq->entry[i];
    }

    if (self->mub > 0)
    {
        // self.T[self.meq:self.m, :self.N_x] = self.A_ub
        for (int i = self->meq; i < self->m; i++)
            for (int j = 0; j < self->N_x; j++)
                self->T->matrix_entry[i][j] = self->A_ub->matrix_entry[i - self->meq][j];
        // b[self.meq:self.m] = self.b_ub
        for (int i = self->meq; i < self->m; i++)
            self->T->matrix_entry[i][self->T->col_size - 1] = self->b_ub->entry[i - self->meq];
    }

    //人工变量与松弛变量系数设置为1
    for (int i = 0; i < self->m; i++)
        for (int j = self->N_x; j < self->N_x + self->m; j++)
            if (i == j - self->N_x)
                self->T->matrix_entry[i][j] = 1.0;

    // T矩阵最后一行表示-c;
    // 第一阶段:
    for (int j = self->N_x; j < self->N_x + self->meq; j++)
        self->T->matrix_entry[self->T->row_size - 1][j] = -1.;
    // for i in range(0, self.meq):
    //     self.T[-1] += self.T[i]
    for (int i = 0; i < self->meq; i++)
        for (int j = 0; j < self->T->col_size; j++)
            self->T->matrix_entry[self->T->row_size - 1][j] += self->T->matrix_entry[i][j];
    self->sign = (int *)malloc(sizeof(int) * (self->m));
    for (int i = 0; i < self->m; i++)
        self->sign[i] = self->N_x + i;
}

void simplex_free(Simplex *self)
{
    vector_free(self->c);
    if (self->b_ub != NULL)
        vector_free(self->b_ub);
    if (self->A_ub != NULL)
        matrix_free(self->A_ub);

    if (self->b_eq != NULL)
        vector_free(self->b_eq);
    if (self->A_eq != NULL)
        matrix_free(self->A_eq);
    if (self->T != NULL)
        matrix_free(self->T);
    if (self->T != NULL)
        free(self->sign);
    free(self);
}

void simplex_solve(Simplex *self)
{
    int num = 0;
    int flag = true;
    int iternum = 0;
    while (flag)
    {

        // 直至所有非基变量检验数小于等于0;
        // 合并多个解的情况,即时非基变量检验数等于0也停止迭代
        // TODO use hack method to optimize speed here
        Vector *temp = vector_alloc(self->T->col_size - 1);
        for (int i = 0; i < temp->size; i++)
            temp->entry[i] = self->T->matrix_entry[self->T->row_size - 1][i];
        double t = vector_max(temp);
        vector_free(temp);
        if (t <= 0)
            flag = false;
        else
        {
            num += 1;
            self->F = simplex_calculate(self);
        }
        if (self->F == 0)
            break;
        iternum++;
        if (iternum >= 5000)
        {
            log_w("iter overflow");
            break;
        }
    }
}

int simplex_calculate(Simplex *self)
{
    //**Hack;
    Vector H;
    H.entry = self->T->matrix_entry[self->T->row_size - 1];
    H.size = self->T->col_size - 1;
    int j_num = vector_argmax(&H);
    //
    Vector *D = vector_alloc(self->m);
    for (int i = 0; i < self->m; i++)
    {
        if (self->T->matrix_entry[i][j_num] == 0.0)
            D->entry[i] = INFINITY;
        else
            D->entry[i] = self->T->matrix_entry[i][self->T->col_size - 1] / self->T->matrix_entry[i][j_num];
    }
    if (vector_max(D) <= 0.0)
        return 0;
    //  i_num = D.index(min([x for x in D if x >= 0]))
    for (int i = 0; i < D->size; i++)
        if (D->entry[i] < 0)
            D->entry[i] = NAN;
    int i_num = vector_argmin(D);
    self->sign[i_num] = j_num;
    double t = self->T->matrix_entry[i_num][j_num];
    //  self.T[i_num] /= t
    for (int j = 0; j < self->T->col_size; j++)
        self->T->matrix_entry[i_num][j] /= t;
    // for i in [x for x in range(0, self.m + 1) if x != i_num]:
    //     self.T[i] -= self.T[i_num] * self.T[i][j_num]
    // TODO fix me DONE
    for (int i = 0; i < self->m + 1; i++)
    {
        if (i != i_num)
        {
            Vector *temp = vector_alloc(self->T->col_size);
            for (int j = 0; j < temp->size; j++)
                temp->entry[j] = self->T->matrix_entry[i_num][j] * self->T->matrix_entry[i][j_num];
            for (int j = 0; j < self->T->col_size; j++)
                self->T->matrix_entry[i][j] -= temp->entry[j];
            vector_free(temp);
        }
        else
            continue;
    }
    vector_free(D);
    return 1;
}

void simplex_change(Simplex *self)
{
    // 人工变量所在列变为0,替换上第二阶段的c
    // self.T[:, self.N_x:self.N_x + self.meq] = 0
    // self.T[-1, 0:self.N_x] = -self.c
    // for i in range(0, self.m):
    //    self.T[-1] -= self.T[i] * self.T[-1][int(self.sign[i])]

    for (int i = 0; i < self->T->row_size; i++)
        for (int j = self->N_x; j < self->N_x + self->meq; j++)
            self->T->matrix_entry[i][j] = 0.0;
    for (int j = 0; j < self->N_x - 1; j++)
    {
        self->T->matrix_entry[self->T->row_size - 1][j] = -self->c->entry[j];
    }
    for (int i = 0; i < self->m; i++)
    {
        Vector *temp = vector_alloc(self->T->col_size);
        for (int j = 0; j < temp->size; j++)
            temp->entry[j] = self->T->matrix_entry[i][j] * self->T->matrix_entry[self->T->row_size - 1][self->sign[i]];
        for (int j = 0; j < self->T->col_size; j++)
            self->T->matrix_entry[self->T->row_size - 1][j] -= temp->entry[j];
        vector_free(temp);
    }
}

int simplex_main(Simplex *self, Vector *xstar)
{
    simplex_inital_value(self);
    if (self->meq > 0)
    {
        log_i("phase 1");
        simplex_solve(self);
        //消除人工变量
        simplex_change(self);
        log_i("pahse 2");
        simplex_solve(self);
    }
    else
    {
        log_i("simple");
        simplex_change(self);
        simplex_solve(self);
    }

    if (self->F == 1)
    {
        log_v("optional solution found");
        int j = 0;
        for (int i = 0; i < self->N_x; i++)
        {
            int flag = false;
            for (int k = 0; k < self->m; k++) // if i in
                if (self->sign[k] == i)
                {
                    flag = true;
                    break;
                }
            if (flag)
            {
                xstar->entry[i] = self->T->matrix_entry[j][self->T->col_size - 1];
                j++;
            }
            else
            {
                xstar->entry[i] = 0.0;
            }
        }
    }
    else
    {
        log_v("infity solution");
    }
    return self->F;
}