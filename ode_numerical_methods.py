import numpy as np
import matplotlib.pyplot as plt

# function for the right-hand side of the differential equation
def f(x, y):
    return np.sin(7 * x) * np.cos(3 * y)**2

# Euler's Method for numerical integration
def euler_method(f, x0, y0, h, xmax):
    x_values = [x0]
    y_values = [y0]
    while x_values[-1] < xmax:
        x = x_values[-1]
        y = y_values[-1]
        y_next = y + h * f(x, y)
        y_values.append(y_next)
        x_values.append(x + h)
    return x_values, y_values

# Modified Euler's Method for numerical integration
def modified_euler_method(f, x0, y0, h, xmax):
    x_values = [x0]
    y_values = [y0]
    while x_values[-1] < xmax:
        x = x_values[-1]
        y = y_values[-1]
        y_predictor = y + h * f(x, y)
        y_corrector = y + h * f(x + h, y_predictor)
        y_next = 0.5 * (y + y_corrector)
        y_values.append(y_next)
        x_values.append(x + h)
    return x_values, y_values

# Runge-Kutta 4th Order Method for numerical integration
def runge_kutta4_method(f, x0, y0, h, xmax):
    x_values = [x0]
    y_values = [y0]
    while x_values[-1] < xmax:
        x = x_values[-1]
        y = y_values[-1]
        k1 = h * f(x, y)
        k2 = h * f(x + 0.5 * h, y + 0.5 * k1)
        k3 = h * f(x + 0.5 * h, y + 0.5 * k2)
        k4 = h * f(x + h, y + k3)
        y_next = y + (1 / 6) * (k1 + 2*k2 + 2*k3 + k4)
        y_values.append(y_next)
        x_values.append(x + h)
    return x_values, y_values

# Adams-Bashforth 4th Order Method for numerical integration
def adams_bashforth_method(f, x0, y0, h, xmax):
    x_values, y_values = runge_kutta4_method(f, x0, y0, h, x0 + 3 * h)
    while x_values[-1] < xmax:
        x = x_values[-1]
        y = y_values[-1]
        y_next = y + h * (55 * f(x, y) - 59 * f(x_values[-2], y_values[-2]) + 37 * f(x_values[-3], y_values[-3]) - 9 * f(x_values[-4], y_values[-4])) / 24
        y_values.append(y_next)
        x_values.append(x + h)
    return x_values, y_values

# Milne's Method for numerical integration
def milne_method(f, x0, y0, h, xmax):
    x_values, y_values = runge_kutta4_method(f, x0, y0, h, x0 + 3 * h)
    while x_values[-1] < xmax:
        x = x_values[-1]
        y = y_values[-1]
        y_predictor = y + 4 * h / 3 * (2 * f(x, y) - f(x_values[-2], y_values[-2]) + 2 * f(x_values[-3], y_values[-3]))
        y_corrector = y + h / 3 * (f(x + h, y_predictor) + 4 * f(x, y) + f(x_values[-2], y_values[-2]))
        y_next = 0.5 * (y_predictor + y_corrector)
        y_values.append(y_next)
        x_values.append(x + h)
    return x_values, y_values

# Initial conditions and step size
x0 = 0
y0 = 0.333
h = 0.01
xmax = 2

# Calculate the numerical solutions
x_euler, y_euler = euler_method(f, x0, y0, h, xmax)
x_modified_euler, y_modified_euler = modified_euler_method(f, x0, y0, h, xmax)
x_runge_kutta, y_runge_kutta = runge_kutta4_method(f, x0, y0, h, xmax)
x_adams_bashforth, y_adams_bashforth = adams_bashforth_method(f, x0, y0, h, xmax)
x_milne, y_milne = milne_method(f, x0, y0, h, xmax)

# Plot all the functions
plt.figure(figsize=(10, 6))
plt.plot(x_euler, y_euler, label='Euler Method')
plt.plot(x_modified_euler, y_modified_euler, label='Modified Euler Method')
plt.plot(x_runge_kutta, y_runge_kutta, label='Runge-Kutta Method')
plt.plot(x_adams_bashforth, y_adams_bashforth, label='Adams-Bashforth Method')
plt.plot(x_milne, y_milne, label='Milne Method')
plt.title("Numerical Solution of y' = sin(7x) * cos^2(3y) with y(0) = 0.333")
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()
