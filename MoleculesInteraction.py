import Dump
import numpy as np
import os


def start(**args):
    AgRadius, AgMass, AlNumOfAtoms, AlRadius, AlMass = \
        args['AgRadius'], args['AgMass'], args['AlNumOfAtoms'], args['AlRadius'], args['AlMass']
    TimeStep, Steps, OutputFrequency, Borders, OutputFileName = \
        args['TimeStep'], args['Steps'], args['OutputFrequency'], args['Borders'], args['OutputFileName']

    dimension = len(Borders)

    AlMass = np.ones(AlNumOfAtoms) * AlMass
    AlRadius = np.ones(AlNumOfAtoms) * AlRadius
    AlVelocity = np.array([[0 for i in range(dimension)] for j in range(AlNumOfAtoms)])
    AlPositions = np.array([[0 for i in range(dimension)] for j in range(AlNumOfAtoms)])

    for i in range(dimension):
        AlPositions[:, i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * AlPositions[:, i]

    if os.path.exists(OutputFileName):
        os.remove(OutputFileName)

    step = 0

    Dump.writeOutput(OutputFileName, AlNumOfAtoms, step, Borders, radius=AlRadius, pos=AlPositions, v=AlVelocity)
