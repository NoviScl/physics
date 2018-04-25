# planetary motion in a plane

import scipy
from scipy.integrate import odeint
from matplotlib import pyplot

def earth(X, t):
    x, y, vx, vy = X
    k = G * M / ((1 + m/M)**2 * (x**2 + y**2)**(3/2))
    return [vx, vy, -k*x, -k*y]
G = 6.67408e-11     # gravitation constant in m^3/(kg*s^2)
m = 5.9722e24       # mass of Earth in kg
M = 1.98855e30      # mass of sun in kg
AU = 1.495978707e11 # astronomical unit in m
day = 86400         # duration of a day in seconds
year = 365.256363004 * day
timeaxis = scipy.linspace(0, 2*year, 100)
x0, y0 = 1.52096e11, 0
vx0, vy0 = 0, 29776.2
#vx0, vy0 = 0, 41000.0
trajectory = odeint(earth, [x0, y0, vx0, vy0], timeaxis)
pyplot.plot(trajectory[:,0], trajectory[:,1], "o-")
pyplot.show()
