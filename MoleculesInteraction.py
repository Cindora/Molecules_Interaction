import Data
import numpy as np


def start(**args):
    AgRadius, AgMass, AlNumOfAtoms, AlRadius, AlMass = \
        args['AgRadius'], args['AgMass'], args['AlNumOfAtoms'], args['AlRadius'], args['AlMass']
    TimeStep, Steps, OutputFrequency, Borders, OutputFileName = \
        args['TimeStep'], args['Steps'], args['OutputFrequency'], args['Borders'], args['OutputFileName']

    AgMass = np.ones(AlNumOfAtoms) * AgMass
    dimension = len(Borders)

    positions = [[0 for i in range(AlNumOfAtoms)] for j in range(dimension)]

    velocity = [[0 for i in range(AlNumOfAtoms)] for j in range(dimension)]

    for i in range(dimension):
        positions[:, i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * positions[:, i]

