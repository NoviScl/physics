# damped mathematical pendulum

import scipy
from scipy.integrate import odeint
from matplotlib import pyplot

g = 9.81    # gravitational acceleration
l = 1.0     # length of pendulum
b = 0.1     # viscous friction coefficient in N*s

def f(X, t):
    phi, omega = X
    alpha = -g/l*scipy.sin(phi) - b*omega
    return [omega, alpha]

x0 = [0, 7]	#initial values(angle and angular velocity)
timeaxis = scipy.linspace(0, 15, 150)

trajectory = odeint(f, x0, timeaxis)

pyplot.title("Mathematical Pendulum (anharmonic)")
pyplot.xlabel("time [s]")
pyplot.ylabel("angle [revolutions]")

pyplot.plot(timeaxis, trajectory[:,0]/(2*scipy.pi))
pyplot.show()


