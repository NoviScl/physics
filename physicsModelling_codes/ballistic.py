# ballistic motion with air resistance

import scipy
from scipy.integrate import odeint
from matplotlib import pyplot

def ballistic(X, t):
    x, y, vx, vy = X
    factor = C / (2*m) * scipy.sqrt((vx-ux)**2 + vy**2)
    return [vx, vy, -factor * (vx-ux), -g - factor * vy]

m = 1.0     # mass in kg
rho = 1.225 # density of air in kg/m^3
A = 0.01    # cross section in m^2
c_w = 0.8   # drag parameter
C = c_w * rho * A
g = 9.81    # gravitational acceleration in m/s^2
ux = -10    # horizontal wind speed in m/s
x0, y0 = 0, 0 # initial position
v = 100.0   # initial speed in m/s

timeaxis = scipy.linspace(0, 15, 100)

degrees = scipy.pi/180 # for converting degrees to radians
for alpha in scipy.linspace(0*degrees, 90*degrees, 10):
    vx, vy = v * scipy.cos(alpha), v * scipy.sin(alpha)
    trajectory = odeint(ballistic, [x0, y0, vx, vy], timeaxis)
    pyplot.plot(trajectory[:,0], trajectory [:,1], "o-")
pyplot.ylim(0)  # show only positive heights
pyplot.show()
