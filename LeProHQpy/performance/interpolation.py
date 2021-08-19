import timeit

import numpy as np
import scipy.interpolate as ip

nx = 10
ny = 10
data = np.random.rand(nx + 1, ny + 1)
data[0, 0] = 0
data[1:, 0] = np.linspace(0, 1, nx)
data[0, 1:] = np.linspace(0, 1, ny)

# print(data)


def recreate_and_call(n=1000):
    def f():
        for j in range(n):
            i = ip.RectBivariateSpline(data[1:, 0], data[0, 1:], data[1:, 1:])
            i(0.5, 0.5)

    return f


def create_once_and_call(n=1000):
    def f():
        i = ip.RectBivariateSpline(data[1:, 0], data[0, 1:], data[1:, 1:])
        for j in range(n):
            i(0.5, 0.5)

    return f


n_it = 100
n_call = 1000

print("recreate", timeit.timeit(recreate_and_call(n_call), number=n_it))
print("create once", timeit.timeit(create_once_and_call(n_call), number=n_it))
