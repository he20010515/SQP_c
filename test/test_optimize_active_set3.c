#include "qp.h"
#include "matrix.h"
#include "vector.h"

int main(int argc, char const *argv[])
{
    sqp_init();
    // Problem
    // f = (x-1)^2 + (y-2)^2 + (z-3)^2 // 这个例子的目标函数
    // s.t. x   == y                   // 等式约束
    //      x+z >= 2                   // 不等式约束

    Matrix *G = matrix_alloc(3, 3); // 申请二次型矩阵的存储空间 G
    double G_array[3][3] = {        // G的元素
                            {
                                2,
                                0,
                                0,
                            },
                            {0, 2, 0},
                            {0, 0, 2}};
    array_2_matrix((double *)G_array, 3, 3, G); //赋值

    Vector *c = vector_alloc(3); // 一次项系数向量
    Vector *b = vector_alloc(2); // 线性约束的右边部分
    c->entry[0] = -2.;           // 赋值
    c->entry[1] = -4;
    c->entry[2] = -6;

    b->entry[0] = 0.;
    b->entry[1] = 2.;

    Matrix *A = matrix_alloc(2, 3); //线性约束矩阵
    double A_array[2][3] = {
        {1, -1, 0},
        {1, 0, 1},
    };
    array_2_matrix((double *)A_array, 2, 3, A);

    LinearConstraints *con = linearconstraints_alloc(3, 2, 1, 1, A, b);
    // 3 问题维度
    // 2 约束数量
    // 1 等式约束数量
    // 1 不等式约束数量
    // A 约束矩阵系数
    // b 约束向量

    Vector *x0 = vector_alloc(3);     //初始点
    Vector *lambda = vector_alloc(2); // lambda

    Vector *xstar = vector_alloc(3); //最优点
    x0->entry[0] = 6;                //初始点赋值 (目前初始点必须是可行解)
    x0->entry[1] = 6;
    x0->entry[2] = 6;

    optimize_qp_active_set(G, c, con, x0, xstar, lambda); //计算
    vector_print(xstar);                                  //打印结果
    return 0;
}
