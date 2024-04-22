import math

# Define functions to be integrated
def f1(x):
    return 1 / math.sqrt(2 * x**2 + 1)

def f2(x):
    return math.log(x + 2) / x

# Trapezoidal Rule for numerical integration
def trapezoidal_rule(func, a, b, n):
    h = (b - a) / n
    integral = 0.5 * (func(a) + func(b))
    for i in range(1, n):
        integral += func(a + i * h)
    integral *= h
    return integral

# Simpson's Rule for numerical integration
def simpsons_rule(func, a, b, n):
    h = (b - a) / n
    integral = func(a) + func(b)
    for i in range(1, n):
        x = a + i * h
        if i % 2 == 0:
            integral += 2 * func(x)
        else:
            integral += 4 * func(x)
    integral *= h / 3
    return integral

# Tolerance for stopping criteria
epsilon = 0.001

# (a) Integrate using Trapezoidal Method
a1, b1 = 0.8, 1.6
n1 = 1
integral_trapezoidal = trapezoidal_rule(f1, a1, b1, n1)
while True:
    n1 *= 2
    new_integral_trapezoidal = trapezoidal_rule(f1, a1, b1, n1)
    if abs(new_integral_trapezoidal - integral_trapezoidal) < epsilon:
        break
    integral_trapezoidal = new_integral_trapezoidal

print("(a) Trapezoidal Method:", integral_trapezoidal)

# (b) Integrate using Simpson's Method
a2, b2 = 1.2, 2
n2 = 1
integral_simpsons = simpsons_rule(f2, a2, b2, n2)
while True:
    n2 *= 2
    new_integral_simpsons = simpsons_rule(f2, a2, b2, n2)
    if abs(new_integral_simpsons - integral_simpsons) < epsilon:
        break
    integral_simpsons = new_integral_simpsons

print("(b) Simpson's Method:", integral_simpsons)
