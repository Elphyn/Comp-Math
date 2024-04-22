import numpy as np

def cubic_spline_interpolation(xi, yi):
    n = len(xi) - 1
    h = [xi[i+1] - xi[i] for i in range(n)]  # Step size between adjacent points
    alpha = [(3/h[i]) * (yi[i+1] - yi[i]) - (3/h[i-1]) * (yi[i] - yi[i-1]) for i in range(1, n)]  # Coefficients for the tridiagonal system
    l = [1] + [2 * (xi[i+1] - xi[i-1]) for i in range(1, n)] + [1]  # Diagonal elements of the tridiagonal matrix
    mu = [0] * (n+1)
    z = [0] * (n+1)
    # Forward elimination
    for i in range(1, n):
        mu[i] = h[i] / (l[i] - mu[i-1] * h[i-1])
        z[i] = (alpha[i-1] - z[i-1] * h[i-1]) / (l[i] - mu[i-1] * h[i-1])
    b = [0] * (n+1)
    c = [0] * (n+1)
    d = [0] * (n+1)
    # Backward substitution
    for j in range(n-1, -1, -1):
        c[j] = z[j] - mu[j] * c[j+1]
        b[j] = (yi[j+1] - yi[j]) / h[j] - h[j] * (c[j+1] + 2 * c[j]) / 3
        d[j] = (c[j+1] - c[j]) / (3 * h[j])  # Coefficients of the cubic polynomial
    return [(yi[i], b[i], c[i], d[i]) for i in range(n)]  # Return spline coefficients

# Given data
xi = [-1.2, -0.9, 0.7, 1.1, 1.7, 2.9]
yi = [3.38688, -1.50579, 16.99677, 25.85121, 28.70127, 0.55419]

# Calculate cubic spline coefficients
spline_coefficients = cubic_spline_interpolation(xi, yi)

# Output spline coefficients
print("Cubic Hermite spline coefficients:")
for i in range(len(spline_coefficients)):
    print(f"For interval [{xi[i]}, {xi[i+1]}]:")
    print("a =", spline_coefficients[i][0])
    print("b =", spline_coefficients[i][1])
    print("c =", spline_coefficients[i][2])
    print("d =", spline_coefficients[i][3])
