import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import csv
import json
import string
import random
import sympy as sp
import numpy as np
import sys
import weightedstats as ws
import pprint
import pandas as pd
from collections import defaultdict
import matplotlib.mlab as mlab
import matplotlib.ticker as mtick

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

These functions are to generate random parameters

"""

def WeightMethod(method):

    options = {
        1 : 2,
        2 : 1
    }

    return options[int(method)]

def WeightChoice (similarity_score):

    options = {
        1: 4,
        2: 2,
        3: 1
    }

    return options[int(similarity_score)]

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

        mode = values[0]

        CI_factor = 10

    else:

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

def ID_Generator(size=10, chars=string.ascii_uppercase + string.digits):
    """
    :param size: the length of the random string
    :param chars: the type of characters in the random string
    :return: a random string to add to our filename as unique ID
    """

    return ''.join(random.choice(chars) for _ in range(size))

def GenerateFileName(data_type, ID):
    """

    :param data_type: the file name in which the data will be saved
    :param ID: the ID of our file
    :return: the filepath in which our file will be saved
    """

    random_part = ID_Generator()
    file_name = ID + "-" + random_part + ".csv"

    file_path = r"C:\\Users\user\PycharmProjects\SimParam\results\\" + data_type + "\\"
    path = file_path + file_name

    return path, file_name

def SaveSamples(param_ID, samples):
    """

    :param param_ID: The ID of the parameter we want to save
    :param samples: the array of samples we generated
    :return: This function will save the file
    """

    file_path, samples_filename = GenerateFileName("Samples", param_ID)

    with open(file_path, 'w', encoding='utf-8-sig') as file:
        file.write(str(param_ID) + ",\n")

        for value in samples:
            file.write(str(value) + ",\n")

    return samples_filename

def SaveMeta(param_ID, mu, sigma, mode, CI_factor, precision):
    """
    :return: This function will save the metadata for 1 parameter
    """

    file_path, meta_filename = GenerateFileName("metadata", param_ID)

    with open(file_path, 'w', encoding='utf-8-sig') as file:
        file.write("mu,sigma,mode,CI factor,precision\n" +
                   str(mu) + "," + str(sigma) + "," + str(mode) + "," + str(CI_factor) + "," + str(precision))

    return meta_filename

def SaveRxn(reaction):
    """

    :param reaction: the reaction object we want to save
    :return: this function will save data for one whole reaction object
    """

    reaction_ID = reaction["ID"]

    file_path, rxn_filename = GenerateFileName("Samples", reaction_ID)

    headers = []
    samples_array = []

    for parameter in reaction["parameters"]:
        headers.append(parameter["ID"])
        samples_array.append(parameter["values"])

    df = pd.DataFrame(samples_array)

    df = df.transpose()
    df.columns = headers

    df.to_csv(path_or_buf=file_path, encoding='utf-8', index=False)

    return rxn_filename

def SaveData(data):
    """
    :param data: the data object we want to save
    :return: This function will save a data object
    """

    rxn_headers = []
    param_headers = []
    samples_array = []

    file_path, data_filename = GenerateFileName("Samples", "total")

    for reaction in data["reactions"]:

        reaction_ID = reaction["ID"]

        for parameter in reaction["parameters"]:
            rxn_headers.append(reaction_ID)
            param_headers.append(parameter["ID"])
            samples_array.append(parameter["values"])

    df = pd.DataFrame(samples_array)
    df = df.transpose()

    df.columns = param_headers

    df.columns = pd.MultiIndex.from_tuples(list(zip(rxn_headers, df.columns)))

    df.to_csv(path_or_buf=file_path, encoding='utf-8', index=False)

    return data_filename

def saveSim (specie_name, specie):

    file_path, samples_filename = GenerateFileName("Sim", specie_name)

    with open(file_path, 'w') as file:

        for value in specie:
            file.write(str(value) + ",\n")

    print("saved at" +file_path)

def saveParam(param_name, values):
    file_path, samples_filename = GenerateFileName("Sim", param_name)

    with open(file_path, 'w') as file:
        for value in values:
            file.write(str(value) + ",\n")

    print("saved at" + file_path)


def SampleParam(values_array, reaction_ID, param_ID, percentage, sample_num):
    """

    this function will coordinate the execution of all other functions and return the samples and the metadata

    """

    values = []
    weights = []

    for value in values_array:
        values.append(float(value[0]))
        weights.append(float(WeightValue(*value[1])))

    mode, CI_factor = CalcModeCI_Factor(values, weights)

    mu, sigma = CalcMuSigma(mode, CI_factor, percentage)

    samples = GenerateSamples(mu, sigma, sample_num)

    return (samples, mu, sigma, mode, CI_factor, 0.00001)

#initializing the global objects


data = '{"reactions":[{"ID":"Glucose Phosphorylation","parameters":[{"ID":"Km_ATP","sampleNum":"10000","percentage":"0.95","values":[["0.03",["1","2","1","3"]],["3.05",["1","2","1","3"]],["2.5",["1","2","1","3"]],["3.1",["1","2","1","3"]]]},{"ID":"Km_gluc","sampleNum":"1000","percentage":"0.95","values":[["6.6",["1","2","1","3"]],["8",["1","2","1","3"]],["8.5",["1","2","1","3"]],["0.3",["1","2","1","3"]],["0.05",["1","2","1","3"]],["1.5",["1","2","1","3"]]]},{"ID":"Km_ADP","sampleNum":"1000","percentage":"0.95","values":[["0.5",["1","2","1","3"]],["0.3",["2","2","1","3"]]]},{"ID":"Km_G6P","sampleNum":"1000","percentage":"0.95","values":[["0.5",["1","2","1","3"]]]},{"ID":"Kcat","sampleNum":"1000","percentage":"0.95","values":[["717",["1","2","1","3"]]]},{"ID":"Keq","sampleNum":"1000","percentage":"0.95","values":[["1310",["1","2","1","3"]]]}]}]}'
#data = '{"reactions":[{"ID":"Glucose Phosphorylation","parameters":[{"ID":"km_glu","sampleNum":"1000","percentage":".95","values":[["0.377",[1,1,1,1]]]},{"ID":"km_ATP","sampleNum":"1000","percentage":".95","values":[["1.84",[1,1,1,1]]]},{"ID":"km_G6P","sampleNum":"1000","percentage":".95","values":[["0.5",[1,1,1,1]]]},{"ID":"km_ADP","sampleNum":"1000","percentage":".95","values":[["0.5",[1,1,1,1]]]},{"ID":"kcat","sampleNum":"1000","percentage":".95","values":[["717",[1,1,1,1]]]},{"ID":"keq","sampleNum":"1000","percentage":".95","values":[["1310",[1,1,1,1]]]}]}]}'
#data = sys.stdin.read()
data = json.loads(data)


for reaction in data["reactions"]:

    reaction_ID = reaction["ID"]

    for parameter in reaction["parameters"]:

        # extract data from incoming object

        param_ID = parameter["ID"]
        percentage = float(parameter["percentage"])
        sample_num = float(parameter["sampleNum"])
        values_array = parameter["values"]

        samples, mu, sigma, mode, CI_factor, precision = SampleParam(values_array, reaction_ID, param_ID, percentage,
                                                                sample_num)

        parameter["mu"] = mu
        parameter["sigma"] = sigma
        parameter["values"] = samples

        #parameter["metadata-filename"] = SaveMeta(param_ID, mu, sigma, mode, CI_factor, precision)

        #parameter["samples-filename"] = SaveSamples(param_ID, samples)

    #reaction["samples-filename"] = SaveRxn(reaction)
#data["samples-filename"] = SaveData(data)

pprint.pprint(data)

print("sampled")

print(data["reactions"][0]["parameters"][0]["values"])

"""
This is to plot parameter values
"""

km_glu = data["reactions"][0]["parameters"][0]["values"]
km_ATP = data["reactions"][0]["parameters"][1]["values"]
km_G6P = data["reactions"][0]["parameters"][2]["values"]
km_ADP = data["reactions"][0]["parameters"][3]["values"]
kcat = data["reactions"][0]["parameters"][4]["values"]
keq = data["reactions"][0]["parameters"][5]["values"]



km_glu_mu = data["reactions"][0]["parameters"][0]["mu"]
km_ATP_mu = data["reactions"][0]["parameters"][1]["mu"]
km_G6P_mu = data["reactions"][0]["parameters"][2]["mu"]
km_ADP_mu = data["reactions"][0]["parameters"][3]["mu"]
kcat_mu = data["reactions"][0]["parameters"][4]["mu"]
keq_mu = data["reactions"][0]["parameters"][5]["mu"]

km_glu_simga = data["reactions"][0]["parameters"][0]["mu"]
km_ATP_simga = data["reactions"][0]["parameters"][1]["mu"]
km_G6P_simga = data["reactions"][0]["parameters"][2]["mu"]
km_ADP_simga = data["reactions"][0]["parameters"][3]["mu"]
kcat_simga = data["reactions"][0]["parameters"][4]["mu"]
keq_simga = data["reactions"][0]["parameters"][5]["mu"]

"""
plt.figure(1)

print("type km glux")
print(type(km_glu))
ax1 = plt.subplot(2,3,1)
plt.title("Km glucose")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(km_glu, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, km_glu_mu, km_glu_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
ax1.locator_params(axis='x', tight=True, nbins=4)
ax1.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)


ax2 = plt.subplot(2,3,2)
plt.title("Km ATP")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(km_ATP, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, km_ATP_mu, km_ATP_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
ax2.locator_params(axis='x', tight=True, nbins=4)
ax2.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)

ax3 = plt.subplot(2,3,3)
plt.title("Km G6P")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(km_G6P, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, km_G6P_mu, km_G6P_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
ax3.locator_params(axis='x', tight=True, nbins=4)
ax3.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)


ax4 = plt.subplot(2,3,4)
plt.title("Km ADP")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(km_ADP, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, km_ADP_mu, km_ADP_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
ax4.locator_params(axis='x', tight=True, nbins=2)
ax4.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)

ax5 = plt.subplot(2,3,5)
plt.title("Keq")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(keq, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, keq_mu, keq_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
#ax5.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
#ax5.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax5.locator_params(axis='x', tight=True, nbins=2)
ax5.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)


ax6 = plt.subplot(2,3,6)
plt.title("Kcat")
plt.ylabel("Frequency")
plt.xlabel("Value/mM")
n, bins, patches = plt.hist(kcat, 50, facecolor='green', alpha=0.75)
#y = mlab.normpdf( bins, kcat_mu, kcat_simga)
#l = plt.plot(bins, y, 'r--', linewidth=1)
#ax6.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
#ax6.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax6.locator_params(axis='x', tight=True, nbins=4)
ax6.locator_params(axis='y', tight=True, nbins=4)
plt.grid(True)

plt.tight_layout()
plt.savefig(r'C:\Users\user\PycharmProjects\SimParam\results\graphs\graphs.png')
#plt.show()

"""


"""
this is to solve the differential equations
"""

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

    # #The Ki is the concentration of inhibitor at which under saturating substrate conditions the reaction rate is half of the maximum reaction rate Vmax
    # #not necessary for random order bi bi
    # #global ki_ATP
    # #global ki_ADP
    #
    # #initial enzyme conc and data
    # hex = 0.02
    # kcat = 717
    #
    # Vf = kcat*hex
    #
    # #setting parameters
    #
    # #ki_ATP = Vf/2
    # #ki_ADP = Vf/2
    #
    # #can you calculate it? no?
    # keq = 1310
    #
    #
    # km_glu = 0.377
    # km_ATP = 1.84
    #
    # #assumed
    # km_G6P = 0.5
    # km_ADP = 0.5


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


    #1 =
    dydt = [-r1, -r1, r1, r1]

    return dydt


#time to run simulation
tspan = np.linspace(0,3600)
#tspan = np.linspace(0,450)

#intitial conditions

y0 = [1, 1, 0, 0]

#defining parameters

def parseData():
    # parameters

    global km_glu
    global km_ATP
    global km_G6P
    global km_ADP

    # enzymatic information
    global hex
    hex = 0.02

    global kcat

    # equilibrium constant
    global keq

    # Vf is forward Vmax, Vr is reverse Vmax
    global Vf

    keq = random.choice(data["reactions"][0]["parameters"][5]["values"])
    kcat = random.choice(data["reactions"][0]["parameters"][4]["values"])
    km_ADP = random.choice(data["reactions"][0]["parameters"][3]["values"])
    km_G6P = random.choice(data["reactions"][0]["parameters"][2]["values"])
    km_ATP = random.choice(data["reactions"][0]["parameters"][1]["values"])
    km_glu = random.choice(data["reactions"][0]["parameters"][0]["values"])
    Vf = kcat * hex

results = []
gluc_60 = []
ATP_60 = []

for _ in range(1000):
    parseData()

    #solving the equation
    solved = odeint(f,y0, tspan)

    # passing/printing results
    #print(solved[:, 0])
    results.append(solved)


plt.figure(2)

for i in range(len(results)):
    # 0 = gluc
    # 1 = ATP
    # 2 = G6P
    # 3 = ADP
    plt.plot(tspan, results[i][:,0], color='b')

    gluc_60.append(results[i][49,0])
    ATP_60.append(results[i][49,1])

# print("gluc_60")
# print(gluc_60)
# print(ATP_60)
#plt.show()


#pprint.pprint(results[0][:,0])
#print(results[0][49,0])
#print(results[0][49,1])


#saveSim("ATP", ATP_60)
saveSim("glucose", gluc_60)

plt.close('all')

print(type(ATP_60))
print(type(gluc_60))
plt.figure(4)

plt.title("ATP concentrations after 7 min for each run")
plt.ylabel("Probability")
plt.xlabel("Value/mM")
weights1 = np.ones_like(ATP_60)/len(ATP_60)
n1, bins1, patches1 = plt.hist(ATP_60, 50, normed=0, facecolor='green', alpha=0.75, weights=weights1)
#n1, bins1, patches1 = plt.hist(ATP_60, 50, facecolor='green', alpha=0.75)


n1 = [float(i)/sum(ATP_60) for i in ATP_60]
print(n1)

# import scipy.stats as st
#
#
# norm = [float(i)/sum(ATP_60) for i in ATP_60]
# # Make a fit to the samples.
# shape, loc, scale = st.lognorm.fit(norm, floc=0)
# # Create an array of length num_bins containing the center of each bin.
# centers = 0.5*(bins1[:-1] + bins1[1:])
# # Compute the CDF at the edges. Then prob, the array of differences,
# # is the probability of a sample being in the corresponding bin.
# cdf = st.lognorm.cdf(bins1, shape, loc=loc, scale=scale)
# prob = np.diff(cdf)
#
# norm = np.asarray(norm)

#plt.plot(centers, norm.size*prob, 'k-', linewidth=2)

plt.grid(True)

plt.savefig(r'C:\Users\user\PycharmProjects\SimParam\results\graphs\ATP_7')
#plt.show()

plt.figure(5)

plt.title("Glucose concentrations after 60min for each run")
plt.ylabel("Probability")
plt.xlabel("Value/mM")
weights2 = np.ones_like(gluc_60)/len(gluc_60)
n2, bins2, patches2 = plt.hist(gluc_60, 50, normed=0, facecolor='green', alpha=0.75, weights=weights2)

# mean_gluc = np.mean(bins2)
# sigma_gluc = np.std(bins2)
# y2 = mlab.normpdf( bins2, mean_gluc, sigma_gluc)
# l2 = plt.plot(bins2, y2, 'r', linewidth=2)

plt.grid(True)

plt.savefig(r'C:\Users\user\PycharmProjects\SimParam\results\graphs\glucose_60')
#plt.show()
