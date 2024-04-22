import numpy as np

def gauss_jordan_elimination(A, B):
    n = len(B)
    # Forward elimination
    for pivot_row in range(n):
        for row in range(pivot_row+1, n):
            multiplier = A[row][pivot_row] / A[pivot_row][pivot_row]
            for col in range(n):
                A[row][col] -= multiplier * A[pivot_row][col]
            B[row] -= multiplier * B[pivot_row]
    # Backward substitution
    for pivot_row in range(n-1, -1, -1):
        for row in range(pivot_row-1, -1, -1):
            multiplier = A[row][pivot_row] / A[pivot_row][pivot_row]
            for col in range(n):
                A[row][col] -= multiplier * A[pivot_row][col]
            B[row] -= multiplier * B[pivot_row]
    # Scaling
    for i in range(n):
        B[i] /= A[i][i]
    return B

def lagrange_interpolation(x_values, y_values, x):
    result = 0
    for i in range(len(x_values)):
        term = y_values[i]
        for j in range(len(x_values)):
            if j != i:
                term *= (x - x_values[j]) / (x_values[i] - x_values[j])
        result += term
    return result

def newton_interpolation(x_values, y_values, x):
    n = len(x_values)
    F = np.zeros((n, n))
    F[:,0] = y_values
    # Divided difference table
    for j in range(1, n):
        for i in range(n-j):
            F[i][j] = (F[i+1][j-1] - F[i][j-1]) / (x_values[i+j] - x_values[i])
    result = F[0][0]
    # Newton's interpolation formula
    for i in range(1, n):
        term = F[0][i]
        for j in range(i):
            term *= (x - x_values[j])
        result += term
    return result

# Given data
x_values = [-5, -3, -1, 1, 3, 5]
y_values = [4, -4, 10, -5, 0, -7]

# Interpolation using Gauss-Jordan elimination method
n = len(x_values)
A = [[x_values[i]**j for j in range(n)] for i in range(n)]
coefficients = gauss_jordan_elimination(A, y_values)
print("Interpolation polynomial coefficients (Gauss-Jordan elimination method):", coefficients)

# User input for interpolation point
user_input = float(input("Enter the interpolation point: "))

# Interpolation using Lagrange's method
lagrange_result = lagrange_interpolation(x_values, y_values, user_input)
print("Interpolated value (Lagrange's method):", lagrange_result)

# Interpolation using Newton's method
newton_result = newton_interpolation(x_values, y_values, user_input)
print("Interpolated value (Newton's method):", newton_result)
