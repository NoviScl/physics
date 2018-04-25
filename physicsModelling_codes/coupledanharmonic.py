# coupled anharmonic oscillator

import scipy
from scipy.integrate import odeint
from matplotlib import pyplot


def coupledpendulums(X, t):
    phi1, phi2, omega1, omega2 = X
    return [omega1, omega2,
            -omega0**2*scipy.sin(phi1) + k*(phi2-phi1), -omega0**2*scipy.sin(phi2) + k*(phi1-phi2)]
omega0 = 1.0
k = 0.001

x0 = [3, 3.01, 0, 0]
#x0 = [3, 3.011, 0, 0]
timeaxis = scipy.linspace(0, 1000, 5000)
trajectory1 = odeint(coupledpendulums, x0, timeaxis)
pyplot.plot(timeaxis, trajectory1[:,0])
pyplot.plot(timeaxis, trajectory1[:,1])
pyplot.show()


