from __future__ import print_function
# import the pysb module and all its methods and functions
from pysb import *
from pysb.integrate import Solver, odesolve
import pylab as pl

#units are in mM

#instantiate the model
Model()

#instantiate species
Monomer('Glucose')
Monomer('G6P')
Monomer('ATP')
Monomer('ADP')

#instantiate parameters
Parameter('Km_Glucose', 0.377)
Parameter('Km_ATP', 1.84)
Parameter('Vmax', 70.75)

#intantiate rules
#Rule('Phosphorylation', None >> Glucose(), Km_Glucose)
Rule('Phosphorylation', Glucose >> G6P , Km_Glucose)
Rule('Adenosine_conversion', ATP >> ADP , Km_ATP)

#Inital conditions
Parameter('Glucose_0', 1)
Parameter('ATP_0', 1)
Initial(Glucose, Glucose_0)
Initial(ATP, ATP_0)

#instanciate Observables

Observable('obsGlucose', Glucose)
Observable('obsG6P', G6P)
Observable('obsATP', ATP)
Observable('obsADP', ADP)


#run simulation
t = pl.linspace(0, 2000)


#solver = Solver(model, t)
#solver.run()

print(solver.y[:, 1])