import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt



def f (y, t):
    """

    :param y: the concentration of the relevant species
    :param t: the in model time
    :return: the solution to the differential equation for glucose, ATP, G6P, ADP
    """

    #parameters

    global km_glu
    global km_ATP
    global km_G6P
    global km_ADP

    #enzymatic information
    global hex
    global kcat

    # equilibrium constant
    global keq

    #Vf is forward Vmax, Vr is reverse Vmax
    global Vf

    #The Ki is the concentration of inhibitor at which under saturating substrate conditions the reaction rate is half of the maximum reaction rate Vmax
    #not necessary for random order bi bi
    #global ki_ATP
    #global ki_ADP

    #initial enzyme conc and data
    hex = 0.02
    kcat = 717

    Vf = kcat*hex

    #setting parameters

    #ki_ATP = Vf/2
    #ki_ADP = Vf/2

    #can you calculate it? no?
    keq = 1310


    km_glu = 0.377
    km_ATP = 1.84

    #assumed
    km_G6P = 0.5
    km_ADP = 0.5


    # Species:

    #species
    #glucose = y(1)
    #ATP = y(2)
    #G6P = y(3)
    #ADP = y(4)

    #rate law (random order bi bi)
    #r1 = ((kcat*hex)*y(1)*y(2)/(km_glu*km_ATP+ y(1) *km_ATP + y(2)*km_glu +y(1)*y(2)))

    #simplifying equation

    alpha_glu = y[0]/km_glu
    alpha_ATP = y[1]/km_ATP
    pi_G6P = y[2]/km_G6P
    pi_ADP = y[3]/km_ADP
    gamma = (y[2]*y[3])/(y[0]*y[1])
    rho = gamma/keq


    r1 = (Vf*(alpha_glu)*(alpha_ATP)*(1-rho)/((1+alpha_glu+pi_G6P)*(1+alpha_glu+pi_ADP)))

    dydt = [-r1, -r1, r1, r1]

    return dydt


#time to run simulation
tspan = np.linspace(0,10)

#intitial conditions

y0 = [1, 1, 0, 0]

#options = odeset('RelTol',1e-9,'AbsTol',1e-14,'NormControl','on')

solve = odeint(f,y0, tspan)

print(solve[:,0])

plt.plot(tspan, solve[:,3], color='b')

plt.show()

# figure(1)
# plot(t,y(:,1), 'b', 'linewidth', 2)
#
# hold on
# plot(t,y(:,2), 'g--')
# plot(t,y(:,3), 'm')
# plot(t,y(:,4), 'c--')
# hold off
#
# title('dunno'), xlabel('Time'),ylabel('welp')
# legend('Glucose', 'ATP', 'G6P', 'ADP')
