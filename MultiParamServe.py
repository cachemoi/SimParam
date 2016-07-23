

import json
import weightedstats as ws
import pprint
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

def weightMethod(method):

    options = {
        1 : 2,
        2 : 1
    }

    return options[method]

def weightChoice (similarity_score):

    options = {
        1: 4,
        2: 2,
        3: 1
    }

    return options[similarity_score]

"""

This method will calculate the total weight of a single value given options

"""


def totalWeight(method, condition, enzyme, organism):

    method_score = weightMethod(method)

    condition_score = weightChoice(condition)
    enzyme_score = weightChoice(enzyme)
    organism_score = weightChoice(organism)

    return condition_score*enzyme_score*organism_score*method_score

"""

This method will find the weighted median of our data

"""

def findWeightedMedian(values_array, total_weights_arrays):

    return ws.weighted_median(values_array, total_weights_arrays)







"""

#initializing the global objects

dataIN = '{"reactions":[{"ID":"lol","param":[{"ID":"","values":[]},{"ID":"","values":[]}]},{"ID":"lol","param":[{"ID":"kkk","values":[["3",["1","2","1","3"]]]},{"ID":"","values":[["3",["1","2","1","3"]],["3",["1","2","1","3"]]]},{"ID":"","values":[]}]}]}'
dataIN = json.loads(dataIN)

dataOUT = DataOUT()


pp = pprint.PrettyPrinter(indent=4)

print(pp.pprint(dataIN))

for key, rxn_array in dataIN.items():

    for rxn in rxn_array:

        rxn_position = rxn_array.index(rxn)

        print(rxn)
        print("new param")

        # initiate the new reaction

        dataOUT.reactions.append(Reaction(rxn["ID"]))

        for param in rxn["param"]:

            print(rxn_position)
            print(param["ID"])

            dataOUT.reactions[rxn_position].parameters.values.append(Parameter(param["ID"]))
"""