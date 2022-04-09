from codecs import BOM_UTF16_BE
from lib2to3.pgen2.token import ATEQUAL
from sre_constants import AT_UNI_BOUNDARY
from matplotlib.image import NonUniformImage
import numpy as np
import time

from scipy.optimize import linprog


class Simplex(object):
    def __init__(self, c, A_ub, b_ub, A_eq, b_eq):  # 初始化函数
        self.c = c  # 系数矩阵
        self.A_ub = A_ub  # 不等式约束矩阵
        self.b_ub = b_ub  # 不等式约束边界
        self.A_eq = A_eq  # 等式约束矩阵
        self.b_eq = b_eq  # 等式约束边界
        self.N_x = 0  # 变量个数
        self.meq = 0  # 等式约束个数
        self.mub = 0  # 不等式约束个数
        self.m = 0  # 约束的总数
        self.T = []  # 将所有数据放入单纯形表T
        self.sign = []  # 记录基变量序号
        self.F = 1  # 判断迭代进行的情况：
        # 1——继续迭代
        # 0——无界解，停止迭代

    def Initial_value(self):
        self.N_x = len(self.c)
        # 按照约束类型求解
        if self.b_eq.any():  # 判断等式约束的有无
            self.meq = len(self.b_eq)
        if self.b_ub.any():  # 判断不等式约束的有无
            self.mub = len(self.b_ub)
        self.m = self.mub + self.meq  # 所有约束的行数
        # 建立第一阶段单纯型表。
        self.T = np.zeros([self.m + 1, self.N_x + self.m + 1], dtype='float')
        b = self.T[:-1, -1]
        # 判断两种约束的有无，
        # T矩阵上半部分为等式约束，下半部分为不等式约束
        if self.meq > 0:
            self.T[:self.meq, :self.N_x] = self.A_eq
            b[:self.meq] = self.b_eq
        if self.mub > 0:
            self.T[self.meq:self.m, :self.N_x] = self.A_ub
            b[self.meq:self.m] = self.b_ub
        # 人工变量与松弛变量系数设为1
        np.fill_diagonal(self.T[:self.m, self.N_x:self.N_x + self.m], 1)
        # T矩阵最后一行表示（-c）
        # 第一阶段
        self.T[-1, self.N_x:self.N_x + self.meq] = -1 * np.ones(self.meq)
        for i in range(0, self.meq):
            self.T[-1] += self.T[i]
        l = [i for i in range(self.N_x, self.N_x + self.m)]  # 记录基变量列标,也可自行输入列标
        self.sign = list(l)

    def solve(self):  # 判断是否到达最优解
        num = 0
        flag = True
        while flag:
            # 直至所有非基变量检验数小于等于0
            # 合并多个解的情况，即使非基变量检验数等于0也停止迭代
            if max(list(self.T[-1][:-1])) <= 0:
                flag = False
            # 迭代，直至某个非基变量检验数等于0
            else:
                num += 1
                self.F = self.calculate()
            # 判断无界解，防止无限迭代
            if self.F == 0:
                break
    # 进行迭代运算

    def calculate(self):
        H = list(self.T[-1, :-1])
        j_num = H.index(max(H))
        H
        D = []
        for i in range(0, self.m):
            if self.T[i][j_num] == 0:
                D.append(float("inf"))
            else:
                D.append(self.T[i][-1] / self.T[i][j_num])
        if max(D) <= 0:  # 判断无界解
            return 0
        # 找出最小比值所在行、列—i_num,j_num
        i_num = D.index(min([x for x in D if x >= 0]))
        self.sign[i_num] = j_num
        t = self.T[i_num][j_num]
        # 换基迭代
        self.T[i_num] /= t
        for i in [x for x in range(0, self.m + 1) if x != i_num]:
            self.T[i] -= self.T[i_num] * self.T[i][j_num]
        return 1

    def change(self):  # 人工变量所在列变为0，替换上第二阶段的c
        self.T[:, self.N_x:self.N_x + self.meq] = 0
        self.T[-1, 0:self.N_x] = -self.c
        for i in range(0, self.m):
            self.T[-1] -= self.T[i] * self.T[-1][int(self.sign[i])]

    def Main(self):  # 主函数
        self.Initial_value()
        if self.meq > 0:  # 采用两阶段法
            # 第一阶段求解
            print("phase 1")
            self.solve()
            print(self.T)
            # 消去人工变量列
            self.change()
            # 第二阶段求解
            print("phase 2")
            self.solve()
        else:
            # 直接进入第二阶段求解，消去人工变量列
            print("simple")
            self.change()
            self.solve()
        if self.F == 1:
            print("Optimal solution:")
            j = 0
            for i in range(0, self.N_x):
                if i not in self.sign:
                    print("x" + str(i + 1) + "=0")
                else:
                    print("x" + str(i+1) + "=", self.T[j][-1])
                    j += 1
            print("Best Value:\n", self.T[-1][-1])
        else:
            print("出现无界解")


if __name__ == '__main__':
    '''
    Min Z=0.4*x1+0.5*x2+0*x3     # Min Z =CX
    st. 0.5*x1+0.5*x2=6             #A_ub*x1<=b_ub
        0.6*x1+0.4*x2-x3=6
        0.3*x1+0.1*x2<=2.7          #A_ub*x1<=b_ub
    '''

    # c = np.array([0.4, 0.5, 0])
    # # 不等式约束<=
    # A_ub = np.array([[0.3, 0.1, 0]])
    # b_ub = np.array([-2.7])
    # # 等式约束
    # A_eq = np.array([[0.5, 0.5, 0],
    #                  [0.6, 0.4, -1]])
    # b_eq = np.array([6, 6])

    # A_ub = np.array([[1, -1, 0],
    #                  [2, 0, -1]])
    # b_ub = np.array([-1, -3])
    # c = np.array([1, -1, 0])
    # A_eq = np.array([1, 0, -1])
    # b_eq = np.array([0.0])

    # x = 4000
    # y = 5*x
    # c = np.random.randint(-10, 10, (y))
    # # 不等式约束<=
    # A_ub = np.array([])
    # b_ub = np.array([])
    # # 等式约束
    # A_eq = np.random.randint(0, 20, (x, y))
    # b_eq = np.random.randint(0, 50, (x))

    A_eq = np.array([[-1, 2, 1, 2, 1, 0, 0, 0, 0],
                     [1, 2, -1, -2, 0, 1, 0, 0, 0],
                     [1, -2, -1, 2, 0, 0, 1, 0, 0],
                     [1, 0, -1, 0, 0, 0, 0, -1, 0],
                     []])

    b_eq = np.array([2, 6, 2, 0, 0])
    A_ub = np.array([])
    b_ub = np.array([])
    c = np.array([-2, -5, 2, 5, 0, 0, 0, 0, 0])
    # print(np.linalg.matrix_rank(A_eq))
    res = linprog(c, None, None, A_eq, b_eq, method='simplex')
    print(res)
    # 输出小数位数
    # time.perf_counter()
    # digit = 2
    # Problem = Simplex(c, A_ub, b_ub, A_eq, b_eq)
    # Problem.Main()
    # print('Running time: %s Seconds', time.perf_counter())
