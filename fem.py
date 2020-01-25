import sys
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt

#parametry równania: funkcje

def a(x):
    return -(x-0.)*(x-0.)-1

def b(x):
    return 1.+4*x

def c(x):
    return -4.

def f(x):
    return 2*(x-0.)*(x-0.)-4*x+3


#parametry równania: stałe

beta = -0.5
gamma = 1.
uR = 0.
n = 10  #liczba elementów skończonych

#definicje parametrów potrzebnych do rozwiazania

h = 1. / n  #podział omegi - przedziału (0, 1) - na n elementów
bignumber = 10000

#funkcje bazowe e

def e(k, x):
    if (k == 0) and (0 <= x) and (x <= h):
        return (h - x) / h
    elif (k == n) and ((n - 1) * h <= x) and (x <= 1):
        return (x / h) - n + 1
    elif (k == n) or (k == 0):
        return 0
    else:  #tutaj wiadomo, że 0<k<n
        x_k = k * h
        if (x <= x_k) and (x_k - h <= x):
            return (x / h) - k + 1
        elif (x <= x_k + h) and (x_k <= x):
            return -(x / h) + k + 1
        else:
            return 0


def ederivative(k, x):
    if (k == 0) and (0 <= x) and (x <= h):
        return -1. / h
    elif (k == n) and ((n - 1) * h <= x) and (x <= 1):
        return 1. / h
    elif (k == n) or (k == 0):
        return 0
    else:  #tutaj wiadomo, że 0<k<n
        x_k = k * h
        if (x <= x_k) and (x_k - h <= x):
            return 1. / h
        elif (x <= x_k + h) and (x_k <= x):
            return -1. / h
        else:
            return 0


def integral1(i, j):
    width = 1. / bignumber
    sum = 0.
    for q in range(bignumber):
        x = q * width
        height = ederivative(j, x) * a(x) * ederivative(i, x)
        area = height * width
        sum += area
    return sum


def integral2(i, j):
    width = 1. / bignumber
    sum = 0.
    for q in range(bignumber):
        x = q * width
        height = e(j, x) * b(x) * ederivative(i, x)
        area = height * width
        sum += area
    return sum


def integral3(i, j):
    width = 1. / bignumber
    sum = 0.
    for q in range(bignumber):
        x = q * width
        height = e(j, x) * c(x) * e(i, x)
        area = height * width
        sum += area
    return sum


def integral4(i):
    width = 1. / bignumber
    sum = 0.
    for q in range(bignumber):
        x = q * width
        height = f(x) * e(i, x)
        area = height * width
        sum += area
    return sum


def Be(i, j):
    return -beta * e(j, 0) * e(i, 0) - integral1(i, j) + integral2(
        i, j) + integral3(i, j)


def le(i):
    return -gamma * e(i, 0) + integral4(i)


# metoda elementów skończonych --- główny program

print("BeMatrix:")
BeMatrix = np.zeros((n + 1, n + 1))
BeMatrix[n, n] = 1
for i in range(n):
    BeMatrix[i, i] = Be(i, i)
    if i > 0:
        BeMatrix[i - 1, i] = Be(i, i - 1)
        BeMatrix[i, i - 1] = Be(i - 1, i)
BeMatrix[n - 1, n] = Be(n, n - 1)
print(BeMatrix)

print("\n", "leVector:")
leVector = np.zeros((n + 1, 1))
for i in range(n):
    leVector[i] = le(i)
leVector[n] = uR
print(leVector)

# rozwiązanie
print("\n", "result:")
result = la.solve(BeMatrix, leVector)
print(result)

# szkic rozwiązania
points = np.linspace(0.0, 1.0, n + 1)
values = np.zeros(n + 1)
for i in range(n + 1):
    values[i] = result[i]
plt.plot(points, result)
plt.show()