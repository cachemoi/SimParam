from sympy import *
from random import *
import sys

def CalcScaleParam_nsolve(mode, CIfact, percentage):

    """
    This function solves the problem with optimization, the problem is that you need a guesstimate.

    We use the function nsolve
    """

    Xmin = mode/CIfact
    Xmax = mode*CIfact

    s = Symbol('s', Real=True)

    eqn = (1/2+1/2*erf((log(Xmax)-(log(mode)+s**2))/(sqrt(2)*s))-(1/2+1/2*erf((log(Xmin)-(log(mode)+s**2))/(sqrt(2)*s)))) - 0.95

    sigma = nsolve(eqn, (0.1,1000000000000000000), solver='bisect')

    mu=log(mode)+sigma**2

    print("sigma=")
    print(sigma)
    print("mu=")
    print(mu.evalf())

def CalcScaleParam_solveset(mode, CIfact, percentage):

    """
    here we use solveset to find all possible solution, it resturns a conditionset with infinity.
    """

    Xmin = mode/CIfact
    Xmax = mode*CIfact

    s = Symbol('s',positive=True)

    eqn = (1/2*erf((log(Xmax)-(log(mode)+s**2))/(sqrt(2)*s))-(1/2*erf((log(Xmin)-(log(mode)+s**2))/(sqrt(2)*s)))) - 0.95

    #returns partial solution ( a set of the possible solutions)

    plot(eqn, (s, 0, 1))

    sigma = solve(eqn, s, domain=S.Reals)

    #mu=log(mode)+sigma**2

    print(sigma)
    #print(mu)

def CalcScaleParam_linsolve (mode, CIfact, percentage):

    #this is an attempt to find a solution with a set. It is a non-linear set of equeation, so linsolve des not work since
    # it only works for a linear set of eqn.

    Xmin = mode / CIfact
    Xmax = mode * CIfact

    s = Symbol('s', Real=True)
    mu = Symbol('mu', Real =True)

    eqn1 = (1/2*erf((log(Xmax)-mu)/sqrt(2)*s)-1/2*erf((log(Xmin)-mu)/sqrt(2)*s)) -percentage
    eqn2 = (exp(mu-s**2)) - mode

    eqns = [eqn1, eqn2]

    plot(eqns)

    ans = linsolve(eqns, s, mu)

    print(ans)

def CalcScaleParam_solve(mode, CIfact, percentage):

    Xmin = mode/CIfact
    Xmax = mode*CIfact

    s = Symbol('s', Real=True)

    eqn = (1/2*erf((log(Xmax)-(log(mode)+s**2))/(sqrt(2)*s))-(1/2*erf((log(Xmin)-(log(mode)+s**2))/(sqrt(2)*s)))) - 0.95

    #returns partial solution ( a set of the possible solutions)

    sigma = solve(eqn, s, set=True)

    #mu=log(mode)+sigma**2

    print(sigma)
    #print(mu)

def CalcScaleParam_nsolveiter(mode, CIfact, percentage):

    """
    this will iterate until a solution is found
    """

    Xmin = mode / CIfact
    Xmax = mode * CIfact

    s = Symbol('s', Real=True)

    eqn = (1 / 2 * erf((log(Xmax) - (log(mode) + s ** 2)) / (sqrt(2) * s)) - (
    1 / 2 * erf((log(Xmin) - (log(mode) + s ** 2)) / (sqrt(2) * s)))) - 0.95

    guess = 0
    should_restart= True

    while should_restart:
        try:
            global sigma
            sigma = nsolve(eqn, guess)
        except:
            guess = guess + 0.1
            print("failed to solve")
        else:
            should_restart = False

    mu = log(mode) + sigma ** 2

    print("sigma=")
    print(sigma)
    print("mu=")
    print(mu.evalf())

def CalcScaleParam_findroot(mode,CIfact, percentage):

    Xmin = mode / CIfact
    Xmax = mode * CIfact

    s = Symbol('s', positive=True)

    eqn = (1 / 2 * erf((log(Xmax) - (log(mode) + s ** 2)) / (sqrt(2) * s)) - (1 / 2 * erf((log(Xmin) - (log(mode) + s ** 2)) / (sqrt(2) * s)))) - 0.95

    eqn = lambdify(s, eqn, 'mpmath')

    sigma= mpmath.findroot(eqn, 0.5, solver="anewton", verbose=True)

    print(sigma)

def CalcScaleParam_bounce(mode,CIfact, percentage, precision, upperbound):

    lowerbound = sys.float_info.min
    testbound=0

    sigma= 0

    Xmin = mode / CIfact
    Xmax = mode * CIfact

    s = Symbol('s')

    eqn = (1 / 2 * erf((log(Xmax) - (log(mode) + s ** 2)) / (sqrt(2) * s)) - (
    1 / 2 * erf((log(Xmin) - (log(mode) + s ** 2)) / (sqrt(2) * s)))) - 0.95

    #plot(eqn, (s, -0.000001, 0.000001), ylim =(-1, 0.5))

    while not ((upperbound -lowerbound) <= precision):
        testbound = lowerbound + (upperbound-lowerbound)/2

        f_lowerbound = eqn.evalf(subs={s: lowerbound})
        f_testbound = eqn.evalf(subs={s: testbound})

        if (f_lowerbound *f_testbound > 0):
            lowerbound=testbound
        else:
            upperbound=testbound
            sigma=testbound

        print("iterated")

    mu = log(mode) + sigma ** 2

    print(mu.evalf())

#CalcScaleParam_solveset(2,2,0.95)
#CalcScaleParam_linsolve(2,2,0.95)
#CalcScaleParam_nsolve(1000000000,100000000,0.99)
#CalcScaleParam_findroot(100,150,0.95)
#CalcScaleParam_nsolveiter(2,2,0.95)
CalcScaleParam_bounce(2,2,0.95, 0.0001, sys.float_info.max)