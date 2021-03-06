

import json

import sympy as sp
import numpy as np
import sys
import weightedstats as ws
from collections import defaultdict

class DataOUT:
    def __init__(self):
        self.reactions = []

class Reaction:
    def __init__(self,ID):
        self.parameters = {}
        self.ID = ID

class Parameter:
    def __init__(self, ID):
        self.ID = ID
        self.samples = []
        self.metadata = {}


"""

These funtions are to weigh the methods and the individual selections

"""

def WeightMethod(method):

    options = {
        1 : 2,
        2 : 1
    }

    return options[method]

def WeightChoice (similarity_score):

    options = {
        1: 4,
        2: 2,
        3: 1
    }

    return options[similarity_score]

def WeightValue(method, condition, enzyme, organism):

    """

    This method will calculate the total weight of a single value given options

    """

    method_score = WeightMethod(method)

    condition_score = WeightChoice(condition)
    enzyme_score = WeightChoice(enzyme)
    organism_score = WeightChoice(organism)

    return condition_score*enzyme_score*organism_score*method_score

def CalcModeCI_Factor (values, weights):

    """

    This function will find the confidence interval factor of our data given an array of values and weights of equal length

    If all values are the same or we only have 1 sample, we use a CI factor of 10 to generate the range of values from which we sample

    """

    if len(values) == 1 or len(set(values)) <= 1:

        print("values were the same")

        mode = values[0]

        CI_factor = 10

        print(mode)
        print(CI_factor)

    else:

        print("ran normal")

        total_weight = sum(weights)
        array_len = len(weights)
        mode = ws.weighted_median(values, weights)

        min_position_weighted = round(0.25*total_weight)
        max_position_weighted = round(0.975*total_weight)

        #translating the weights and the values array to percentages and using this to find out where the min and max positiions are

        min_position= round(((min_position_weighted*100)/total_weight)*array_len/100) -1
        max_position= round(((max_position_weighted*100)/total_weight)*array_len/100) -1

        sort = sorted(values)

        min_value = sort[min_position]
        max_value = sort[max_position]

        if min_value == max_value:

            print("weighted values were the same")

            mode = min_value
            CI_factor = 10

        else:
            CI_min = mode/min_value
            CI_max = max_value/mode

            CI_factor = (CI_min + CI_max)/2

    return (mode, CI_factor)

def CalcMuSigma(mode,CIfact, percentage):

    """


    This will calculate the value of sigma and mu of a theoretical log normal distribution with an arbitrary amount of
    precision given a Confidence interval factor, a mode and a percentage of values that will be contained within the defined boundaries

    The equation to find mu and sigma is a system of 2 equation.
        - The CDF
        - The mode

    The system has been rewritten so that we can find the sigma first, and then use this to determine the mu.

    """
    upperbound = sys.float_info.max
    lowerbound = sys.float_info.min
    precision = 0.00001

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

    mu = sp.log(mode) + sigma ** 2

    mu = float(mu.evalf())

    return (mu, sigma)

def GenerateSamples(mu, sigma, sample_size):

    """

    :param mu:
     the mean of a population

    :param sigma:
    measure of variation and dispesion of a set of data value

    :sample_size:
    The desired number of samples

    :return:
    This function will return the desired number of random samples taken from the lognormal distribution with the
    specified mu and sigma.

    It also outputs a graph for now (debugging)

    """

    samples = np.random.lognormal(mu, sigma, sample_size)
    samples = samples.tolist()

    return samples

def getSamples(values_array,  reaction_ID, param_ID, percentage, sample_num):

    values = []
    weights = []

    for value in values_array:

        values.append(float(value[0]))
        weights.append(float(WeightValue(*value[1])))

    mode, CI_factor = CalcModeCI_Factor(values,weights)

    mu, sigma = CalcMuSigma(mode, CI_factor, percentage)

    samples = GenerateSamples(mu, sigma, sample_num)

    return samples

#initializing the global objects


data = sys.stdin.read()
data = json.loads(data)


for reaction in data["reactions"]:

    reaction_ID = reaction["ID"]

    for parameter in reaction["parameters"]:

        #extract data from incoming object

        param_ID = parameter["ID"]
        percentage = float(parameter["percentage"])
        sample_num = float(parameter["sampleNum"])
        values_array = parameter["value"]

        samples = getSamples(values_array, reaction_ID, param_ID, percentage, sample_num)

        parameter["value"] = samples

print(json.dumps(data))