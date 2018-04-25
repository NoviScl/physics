# Demonstration of numerical integration algorithms
# using harmonic oscillator as an example

import scipy
from matplotlib import pyplot

def EulerStep(f, x0, t0, t1):
    "Single Euler integration step"
    h = t1 - t0
    k1 = f(x0, t0)
    return h * k1

def HeunStep(f, x0, t0, t1):
    "Single Heun integration step"
    h = t1 - t0
    k1 = f(x0, t0)
    k2 = f(x0 + h*k1, t0+h)
    return h * (k1 + k2)/2

def RK4Step(f, x0, t0, t1):
    "Single 4-th order Runge-Kutta integration step"
    h = t1 - t0
    k1 = f(x0, t0)
    k2 = f(x0 + h/2*k1, t0+h/2)
    k3 = f(x0 + h/2*k2, t0+h/2)
    k4 = f(x0 + h*k3, t0+h)
    return h * (k1 + 2*k2 + 2*k3 + k4)/6

def solveODE(integrator, f, x0, timeaxis):
    "Numerical integration of ODE"
    # convert the list of initial values (e.g. x, v) to a vector
    x = scipy.array(x0)
    # list holding integrated data points, starting with the initial value
    sums = [x]
    # apply integration step method and record the results
    for i in range(1, len(timeaxis)):
        t0 = timeaxis[i-1]
        t1 = timeaxis[i]
        x = x + integrator(f, x, t0, t1)
        sums.append(x)
    return scipy.array(sums)

# Physical system: spring-mass system
def f(X, t):
    # unpack the 2-dimensional vector X = (x, v) into its components
    x, v = X
    a = -k/m * x # Hooke's law: F = - k * x
    return scipy.array([v, a]) # repack result to a 2-dimensional vector

# Set up the parameters defining the system
k = 10      # spring constant in N/m
m = 1.0     # mass in kg
x0 = 1.0    # initial elongation in m
v0 = 0.5    # initial velocity in m/s
X0 = [x0, v0]   # initial state of the system at the start time
startTime = 0 # in s
stopTime = 5 # in s
steps = 50

# generate a sequence of times from startTime ... stopTime in 'steps' steps
timeaxis = scipy.linspace(startTime, stopTime, steps)

# Prepare the graph
pyplot.title("Effect of Integration Method (all using %d steps)"%steps)
pyplot.ylabel("elongation")
pyplot.xlabel("time")

# Integrate the ODE using different methods and plot the results
trajectory = solveODE(EulerStep, f, X0, timeaxis)
pyplot.plot(timeaxis, trajectory[:, 0], color = "blue")
pyplot.text(0.25, 4, "Euler", color="blue")

trajectory = solveODE(HeunStep, f, X0, timeaxis)
pyplot.plot(timeaxis, trajectory[:, 0], color = "green")
pyplot.text(0.25, 5, "Heun", color="green")

trajectory = solveODE(RK4Step, f, X0, timeaxis)
pyplot.plot(timeaxis, trajectory[:, 0], color = "red")
pyplot.text(0.25, 6, "RK4", color="red")

# Plot the analytical (exact) result for comparison
omega = scipy.sqrt(k/m)
trueTrajectory = X0[0] * scipy.cos(omega * timeaxis) + X0[1]/omega * scipy.sin(omega * timeaxis)
pyplot.plot(timeaxis, trueTrajectory, "black")
pyplot.text(0.25, 7, "exact solution", color="black")

pyplot.show()

