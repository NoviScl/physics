# Demonstration of effect of step size on numerical integration
# using harmonic oscillator as an example

import scipy
from matplotlib import pyplot

def EulerStep(f, x, t, h):
    "Single Euler integration step"
    return x + h * f(x, t)

def integrate(stepMethod, f, x0, timeaxis):
    "Numerical integration of ODE using chosen method"
    x = scipy.array(x0)
    result = [x]
    t = timeaxis[0]
    for tNext in timeaxis[1:]:
        h = tNext - t
        x = stepMethod(f, x, t, h)
        result.append(x)
        t = tNext
    return scipy.array(result)

# Physical system: spring-mass system
k = 10      # spring constant in N/m
m = 1.0     # mass in kg

def f(X, t):
    x, v = X
    a = -k/m * x
    return scipy.array([v, a])

initialElongation = 1.0
initialVelocity = 5.0

x0 = [initialElongation, initialVelocity]

startTime = 0
stopTime = 5

pyplot.title("Effect of Integration Step Size")
pyplot.ylabel("elongation")
pyplot.xlabel("time")

textpos = 5
colors = ["blue", "green", "red", "black"]

for steps in 50, 250, 1000:
    timeaxis = scipy.linspace(startTime, stopTime, steps)
    trajectory = integrate(EulerStep, f, x0,timeaxis)
    color = colors.pop(0)
    pyplot.plot(timeaxis, trajectory[:,0], color=color)
    pyplot.text(0.25, textpos, "%d steps"%steps, color=color)
    textpos += 1

# analytical solution for comparison
omega = scipy.sqrt(k/m)
trueTrajectory = x0[0] * scipy.cos(omega * timeaxis) \
                 + x0[1]/omega * scipy.sin(omega * timeaxis)

color=colors.pop(0)
pyplot.plot(timeaxis, trueTrajectory, color=color)
pyplot.text(0.25, textpos, "exact solution")
pyplot.show()

