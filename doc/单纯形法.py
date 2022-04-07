from fractions import Fraction as f

import numpy as np


def getinput():
    global m, n  # 这两个变量其他函数里也需要调用
    # string = input('''
    # 输入初始单纯形表形如
    # 例一：3 2 1 0 18;-1 4 0 1 8;-2 1 0 0 0
    # 例二：2 1 0 1 0 0 8;-4 -2 3 0 1 0 14;1 -2 1 0 0 1 18;6 -3 1 0 0 0 0
    # 例三：8 2 4 1 0 0 1;2 6 6 0 1 0 1;6 4 4 0 0 1 1;1 1 1 0 0 0 0
    # 前m行表示m个约束的增广矩阵，最后一行表示检验数[价值向量]
    # 输入：''')
    string = "1 0 -1 0 0 1 0 0 0;0 1 0 -1 0 0 1 0 0;1 1 1 0 1 0 0 1 1;0 0 0 0 0 -1 -1 -1 0"
    a = [list(map(eval, row.split())) for row in string.split(';')]
    matrix = np.array(a)
    m, n = matrix.shape
    n -= 1
    m -= 1
    print('\n\n输入的目标函数为')
    x = [f'{matrix[-1,j]}*x_{j+1}' for j in range(n)]
    print('max z = '+' + '.join(x))
    print('\n\n输入的方程为')
    for i in range(m):
        x = [f'{matrix[i,j]}*x_{j+1}' for j in range(n)]
        print(' + '.join(x), f'={matrix[i,-1]}')
    print(f'\n\n有{m}个约束条件，{n}个决策变量')
    return a


def judge(matrix):
    if max(matrix[-1][:-1]) <= 0:  # 最后一行除了b列的所有检验数
        flag = False
    else:
        flag = True
    return flag


def pr(matrix, vect):  # 输出单纯形表
    print('*'*20)
    print(' ', end='\t')
    for i in range(n):
        print('X_{}'.format(i+1), end='\t')
    print('b')
    for i in range(m+1):
        if i <= m-1:
            print('x_{}'.format(vect[i]), end='\t')
        elif i == m:
            print('r_1', end='\t')
        for j in matrix[i]:
            print(f(str(j)).limit_denominator(), end='\t')  # 输出分数形式
        print(end='\n')


def trans(a, matrix, vect):  # 转轴
    maxi = max(matrix[-1][:-1])
    index = a[-1].index(maxi)  # 入基变量的足标

    l = {}
    for i in a[:-1]:  # 遍历矩阵的每一行
        if i[index] > 0:  # 如果这一行的第index个元素大于0了
            l[i[-1]/i[index]] = a.index(i)  # 记录下来当前行的坐标
    pivot = l[min(l)]  # 出基变量的足标

    matrix[pivot] = matrix[pivot]/matrix[pivot][index]

    for i in range(len(a)):
        if i != pivot:
            matrix[i] = matrix[i] - matrix[i][index]*matrix[pivot]

    vect[pivot] = index+1  # 基变量足标变化
    a = [list(i) for i in matrix]  # 把原来的列表也同时变换掉，为了方便索引

    return a, matrix, vect


def print_solution(matrix, vect):
    print('*'*20)
    for i in range(1, n+1):
        if i in vect:
            print('x{}*={}'.format(i,
                  f(str(matrix[vect.index(i)][-1])).limit_denominator()), end='，')
        else:
            print('x{}*={}'.format(i, 0), end='，')
    print('\nz*={}'.format(f(str(-matrix[-1][-1])).limit_denominator()))


def main():
    a = getinput()
    matrix = np.array(a, dtype=np.float64)  # array可以方便地进行整行的操作，而列表可以方便索引
    # vect = [int(input('输入基变量足标')) for i in range(m)]
    vect = [4, 5, 6]
    pr(matrix, vect)
    while judge(matrix):
        a, matrix, vect = trans(a, matrix, vect)
        # pr(matrix, vect)
        print(vect)
    print_solution(matrix, vect)


if __name__ == '__main__':
    # main()

    A = np.array([
        [11, 12, 13, 14],
        [21, 22, 23, 24],
        [31, 32, 33, 34],
        [41, 42, 43, 44],
    ])
    idx = np.array([1, 2, 4, 3], dtype=int) - 1
    print(A[idx, :][:, idx])
