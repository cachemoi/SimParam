import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def f(y,t):
    dydt = y*np.exp(-t)
    return dydt

#initial conditions on y
y0 = [1]

y = odeint(f, y0, [0,0.2,0.5] )

# instantiate the parameters
#vmax = 70.75
#kma = 0.3
#kmb = 1.84

#concentrations
#A = 1
#B = 1

#calculate reaction rate (bi-bi rate law)
#v = (vmax*A*B)/(kma*kmb+ A*kmb + B*kma + A*B)

#print(v)