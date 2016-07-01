
from sympy import *


def CalcScaleParamOPT(mode, percentage, Xmin, Xmax):

    s = Symbol('s', Real=True)

    eqn = (1/2+1/2*erf((log(Xmax)-(log(mode)+s**2))/(sqrt(2)*s))-(1/2+1/2*erf((log(Xmin)-(log(mode)+s**2))/(sqrt(2)*s)))) - 0.95

    sigma = solveset(eqn, s)

    #mu=log(mode)+sigma**2

    print(sigma)
    #print(mu.evalf())

def CalcScaleParamSET (mode, percentage, Xmin, Xmax):

    s = Symbol('s', Real=True)
    mu = Symbol('mu', Real =True)

    eqn1 = (1/2+1/2*erf((log(Xmax)-mu)/sqrt(2)*s)-1/2-1/2*erf((log(Xmin)-mu)/sqrt(2)*s)) -percentage
    eqn2 = (exp(mu-s**2)) - mode

    eqns = [eqn1, eqn2]

    ans = linsolve(eqns, (s, mu))

    print(ans)

def CalcScaleParam_solveset(mode, percentage, Xmin, Xmax):

    s = Symbol('s', positive=True)

    eqn = (1 / 2 * erf((log(Xmax) - (log(mode) + s ** 2)) / (sqrt(2) * s)) - (1 / 2 * erf((log(Xmin) - (log(mode) + s ** 2)) / (sqrt(2) * s)))) - 0.95

    plot(eqn, (s,5059707,65059707))

    eqn = lambdify(s, eqn, 'mpmath')

    #sigma = solveset(eqn, s, domain=S.Reals)

    sigma= mpmath.findroot(eqn)

    print(sigma)


#CalcScaleParam_solveset(4, 0.95, 5, 6)
CalcScaleParam_solveset(2, 0.95, 1, 4)

#CalcScaleParamOPT(2, 0.95, 1, 4)

#CalcScaleParamSET(2, 0.95, 1, 4)