import numpy as np
import scipy


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


def SQP():
    scipy.
    pass


if __name__ == "__main__":
    pass
