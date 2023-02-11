import Dump
import numpy as np
import os


def calcForces(Masses, Velocities, TimeStep):
    numOfAtoms, dimension = Velocities.shape

    someForce = np.random.randn(numOfAtoms, dimension) * np.sqrt(2.0 * Masses / TimeStep)[np.newaxis].T
    forces = (Velocities * Masses[np.newaxis].T) + someForce

    return forces


def applyForces(Positions, Velocities, forces, Masses, TimeStep):
    Positions += Velocities * TimeStep
    Velocities += forces * TimeStep / Masses[np.newaxis].T



def start(**args):
    """ Считывание данных из кортежа """
    AgRadius, AgMass, AlNumOfAtoms, AlRadius, AlMass = \
        args['AgRadius'], args['AgMass'], args['AlNumOfAtoms'], args['AlRadius'], args['AlMass']
    TimeStep, Steps, OutputFrequency, Borders, OutputFileName = \
        args['TimeStep'], args['Steps'], args['OutputFrequency'], args['Borders'], args['OutputFileName']

    dimension = len(Borders)  # Вычислние размерности

    """ Преобразование значений в массивы для удобного использования """

    AgPosition = np.array([0.0 for i in range(dimension)])
    AgVelocity = np.array([0.0 for i in range(dimension)])

    AlMass = np.ones(AlNumOfAtoms) * AlMass
    AlRadius = np.ones(AlNumOfAtoms) * AlRadius
    AlVelocity = np.array([[0.0 for i in range(dimension)] for j in range(AlNumOfAtoms)])
    AlPositions = np.array([[0.0 for i in range(dimension)] for j in range(AlNumOfAtoms)])

    AgPosition = np.array([1/2, 0, 0])

    if AlNumOfAtoms % 2 == 0:
        coeff = 1 / (AlNumOfAtoms // 2)
    else:
        coeff = 1 / ((AlNumOfAtoms+1) // 2)

    for i in range(0, AlNumOfAtoms):
        AlPositions[i][0], AlPositions[i][1] = 1 / 4 + i % 2 * 1 / 2, i // 2 * coeff + coeff

    for i in range(dimension):  # Подгонка значений относительно границ области
        AlPositions[:, i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * AlPositions[:, i]
        AgPosition[i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * AgPosition[i]

    Positions = np.append([AgPosition], AlPositions, axis=0)
    Velocities = np.append([AgVelocity], AlVelocity, axis=0)
    Masses = np.append([AgMass], AlMass)
    Radius = np.append([AgRadius], AlRadius)

    if os.path.exists(OutputFileName):
        os.remove(OutputFileName)

    step = 0

    while step <= Steps:
        step += 1

        forces = calcForces(Masses, Velocities, TimeStep)
        applyForces(Positions, Velocities, forces, Masses, TimeStep)

        Dump.writeOutput(OutputFileName, AlNumOfAtoms + 1, step, Borders,
                         radius=Radius, pos=Positions, v=Velocities)  # Вывод данных в .dump файл
