#This algorithm is written to process the heauristics data
#in Option Object:
#  v = values
#  eq=equal
#  plaus = plausible
# all statements of similarity are translated into 3 integers:
# 1 = same
# 2 = similar
# 3 = different

from collections import defaultdict
import weightedstats as ws

class Data:
    def __init__(self):
        self.value_options = defaultdict(list)
        self.measured_v = []
        self.weights = []
        self.v_eq_plaus = False

    def reset(self):
        self.value_options = defaultdict(list)
        self.measured_v = []
        self.weights = []
        self.v_eq_plaus = False


class Options:
    def __init__ (self):
        self.method = ""
        self.organism = ""
        self.enzyme = ""
        self.species = ""
        self.condition = ""
    def reset(self):
        self.method = ""
        self.organism = ""
        self.enzyme = ""
        self.species = ""
        self.condition = ""

class Weighting :
    def __init__(self):
        self.method = 0
        self.organism = 0
        self.enzyme = 0
        self.conditions = 0
        self.total = 0
    def reset (self):
        self.method = 0
        self.organism = 0
        self.enzyme = 0
        self.conditions = 0
        self.total = 0

options = Options()
data = Data()
weight = Weighting()

def weightMethod(method):

    def in_vivo():
        #if same as modelling system???
        weight.method += 2
    def in_vitro():
        weight.method -= 1

    weight_options = {
        1 : in_vivo,
        2 : in_vitro
    }

    weight_options[method]()


def weightOrganism (organism):

    def same_organism():
        weight.organism = 4
    def phylo_related():
        weight.organism = 2
    def unrelated():
        weight.organism = 1

    organism_options = {
        1 : same_organism,
        2 : phylo_related,
        3 : unrelated
    }

    organism_options[organism]()

#it can be anything, not only enzyme, do I need to deal with more than one case or is this a single value?

def weightEnzyme(enzyme):
    def same():
        weight.enzyme = 4
    def related():
        weight.enzyme = 2
    def unrelated():
        weight.enzyme = 1

    enzyme_options = {
        1 : same,
        2 : related,
        3 : unrelated
    }

    enzyme_options[enzyme]()

def weightConditions(conditions):
    def same_pH_temp():
        weight.conditions = 4
    def similar_pH_temp():
        weight.conditions = 2
    def unrelated():
        weight.conditions = 1

    cond_options ={
        1 : same_pH_temp,
        2 : similar_pH_temp,
        3 : unrelated
    }
    cond_options[conditions]()

def multiplyWeights():

    weight.total = weight.conditions*weight.enzyme*weight.organism*weight.method

def weightValues(value, method, condition, enzyme, organism):

    weightMethod(method)
    weightConditions(condition)
    weightEnzyme(enzyme)
    weightOrganism(organism)

    multiplyWeights()

def main():

    #this is an array with possible parameter values and weight used as example

    data.value_options = {
        2 : [1,2,3,1],
        10 : [1,3,1,2],
        3 : [2,3,1,2],
        5 : [2,3,2,1],
        11 : [1,2,3,2]
    }

    for value in data.value_options:

        data.measured_v.append(value)

        method = data.value_options[value][0]
        condition = data.value_options[value][1]
        enzyme = data.value_options[value][2]
        organism = data.value_options[value][3]


        weightValues(value , method , condition , enzyme , organism)

        data.weights.append(weight.total)

    print(data.weights)
    print(data.measured_v)

    print(ws.median(data.measured_v))


    print(ws.weighted_median(data.measured_v, data.weights))

main()

#most weights can be divided in 3 statements: same, similar, different
#why are we looking for the weighted median? to define the distribution? why not mean?
#species mean chem species, of which enzyme are made of
#make one general specie opriton