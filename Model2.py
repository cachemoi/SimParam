import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


def f(y,t, km_gluc, km_ATP):

    gluc_0, ATP_0 = y
    km_gluc, km_ATP = km_gluc, km_ATP

    dydt = [gluc_0, -km_gluc*ATP_0 - km_ATP*np.sin(gluc_0)]

    return dydt

#initial conditions on y
y0 = 1, 1, 0, 0

#not needed?
gluc_0 = 1
ATP_0 = 1

#parameters
km_gluc = 0.377
km_ATP = 1.84
vmax = 70.75


time = np.linspace(0,60,120)

y = odeint(f, y0, time, args=(km_gluc, km_ATP))

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