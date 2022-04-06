from telnetlib import X3PAD
import numpy as np
import scipy

NUMERICAL_DIFF_STEP = 0.001


def targetfunction(X: np.ndarray):
    x = X[0]
    y = X[1]
    z = X[2]
    return np.exp((x-1)**2 + (y-2)**2 + (z-3)**3)


def con(X: np.ndarray):
    x = X[0]
    y = X[1]
    z = X[2]
    return np.array([x**2 + y**2 + z**2 - 9, x-1, y-1, z-1])


def grad(f, X: np.ndarray):
    g = np.ones_like(X)
    for i in range(X.shape[0]):
        temp = X
        temp[i] += NUMERICAL_DIFF_STEP
        fah = f(temp)
        temp = X
        temp[i] -= NUMERICAL_DIFF_STEP
        fsh = f(temp)
        g[i] = (fah-fsh)/(2*NUMERICAL_DIFF_STEP)
    return g


def SQP(fun, con, x0, lambda0):
    rho = 0.5
    aita = 0.25
    tao = 0.5
    pass


if __name__ == "__main__":
    A = set([(1, 2), (1, 3), (2, 1), (2, 2), (3, 1)])

    X = list(set([x for x, y in A]))
    Y = [sum([b for (a, b) in A if a == x])
         for x in X]
    print(X, Y)
    pass
