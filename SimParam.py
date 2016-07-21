import sympy as sp
import sys
import numpy as np
import matplotlib.pyplot as plt

def CalcScaleParam_bounce(mode,CIfact, percentage, precision):
    """

    This will calculate the value of sigma and mu of a theoretical log normal distribution with an arbitrary amount of
    precision given a Confidence interval factor and a mode.

    The equation to find mu and sigma is a system of 2 equation.
        - The CDF
        - The mode

    The system has been rewritten so that we can find the sigma first, and then use this to determine the mu.

    """

    upperbound = sys.float_info.max
    lowerbound = sys.float_info.min

    Xmin = mode / CIfact
    Xmax = mode * CIfact

    s = sp.Symbol('s')

    eqn = (1 / 2 * sp.erf((sp.log(Xmax) - (sp.log(mode) + s ** 2)) / (sp.sqrt(2) * s)) - (
    1 / 2 * sp.erf((sp.log(Xmin) - (sp.log(mode) + s ** 2)) / (sp.sqrt(2) * s)))) - percentage

    while not ((upperbound -lowerbound) <= precision):
        testbound = lowerbound + (upperbound-lowerbound)/2

        f_lowerbound = eqn.evalf(subs={s: lowerbound})
        f_testbound = eqn.evalf(subs={s: testbound})

        if (f_lowerbound *f_testbound > 0):
            lowerbound=testbound
        else:
            upperbound=testbound
            sigma= float(testbound)
        print('iterating')

    mu = sp.log(mode) + sigma ** 2

    mu = float(mu.evalf())

    return (mu, sigma)

def MyLogRand(mu, sigma, sample_size):

    """

    :param mu:
     μ, the mean of a population

    :param sigma:
    σ, measure of variation and dispesion of a set of data value

    :sample_size:
    The desired number of samples

    :return:
    This function will return the desired number of random samples taken from the lognormal distribution with the
    specified mu and sigma.

    It also outputs a graph for now (debugging)

    """

    sample = np.random.lognormal(mu, sigma, sample_size)

    count, bins, ignored = plt.hist(sample, 100, normed=True, align='mid')
    x = np.linspace(min(bins), max(bins), 10000)
    pdf = (np.exp(-(np.log(x) - mu) ** 2 / (2 * sigma ** 2))
           / (x * sigma * np.sqrt(2 * np.pi)))

    plt.plot(x, pdf, linewidth=2, color='r')
    plt.xlim(xmax=20)
    # plt.axis('tight')
    plt.show()

if __name__ == '__main__':

    mu, sigma = CalcScaleParam_bounce(2,2,0.95, 0.0001)

    print(mu, sigma)

    MyLogRand(mu, sigma, 100000)