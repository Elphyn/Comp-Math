import numpy as np

#functions f(x) for equations
def f1(x):
    return 4*x + 2**x + 6

def f2(x):
    return -2*x**3 - 2*x**2 + 2*x + 5

# Bisection method for root refinement
def bisection_method(f, a, b, tol):
    iterations = 0
    while (b - a) / 2 > tol:
        c = (a + b) / 2
        if f(c) == 0:
            return c, iterations
        elif np.sign(f(c)) == np.sign(f(a)):
            a = c
        else:
            b = c
        iterations += 1
    return (a + b) / 2, iterations

# Hybrid method combining chord and tangent methods for root refinement
def hybrid_method(f, f_prime, a, b, tol):
    iterations = 0
    while abs(b - a) > tol:
        x = (a*f(b) - b*f(a)) / (f(b) - f(a))
        if f(x) == 0:
            return x, iterations
        elif np.sign(f(a)) * np.sign(f(x)) < 0:
            b = x
        else:
            a = x
        iterations += 1
    return (a + b) / 2, iterations

# Parameters for root refinement
tolerance = 1e-6
a1, b1 = -10, 10
a2, b2 = -10, 10

# Finding intervals where functions change sign
roots_intervals1 = []
roots_intervals2 = []

for i in range(a1, b1):
    if f1(i) * f1(i+1) < 0:
        roots_intervals1.append((i, i+1))

for i in range(a2, b2):
    if f2(i) * f2(i+1) < 0:
        roots_intervals2.append((i, i+1))

# Refining roots with the given accuracy and counting iterations
for interval in roots_intervals1:
    root_bisection, iterations_bisection = bisection_method(f1, interval[0], interval[1], tolerance)
    root_hybrid, iterations_hybrid = hybrid_method(f1, lambda x: 4 + np.log(2) * 2**x, interval[0], interval[1], tolerance)
    print("Root in interval {}: Bisection Method: {}, Iterations: {}, Hybrid Method: {}, Iterations: {}".format(interval, root_bisection, iterations_bisection, root_hybrid, iterations_hybrid))

for interval in roots_intervals2:
    root_bisection, iterations_bisection = bisection_method(f2, interval[0], interval[1], tolerance)
    root_hybrid, iterations_hybrid = hybrid_method(f2, lambda x: -6 * x**2 - 4 * x + 2, interval[0], interval[1], tolerance)
    print("Root in interval {}: Bisection Method: {}, Iterations: {}, Hybrid Method: {}, Iterations: {}".format(interval, root_bisection, iterations_bisection, root_hybrid, iterations_hybrid))
