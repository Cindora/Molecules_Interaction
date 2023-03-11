import Dump
import numpy as np
import math
import os

AgEpsilon = 0.00801  
AgSigma = 3.54
AlEpsilon = 0.03917
AlSigma = 2.83
AlEps = 0.27
AlAlpha = 1.16

def calcForces(Velocities, Positions):
    numOfAtoms, dimension = Velocities.shape

    forces = np.array([[0.0 for i in range(dimension)] for j in range(numOfAtoms)])

    '''Вычисление силы из потенциала Леннарда_Джонса'''
    for i in range(1, numOfAtoms):
        R = math.sqrt((Positions[0][0] - Positions[i][0]) ** 2 +
                      (Positions[0][1] - Positions[i][1]) ** 2 +
                      (Positions[0][2] - Positions[i][2]) ** 2)
        Force: float
        Force = (48 * AgEpsilon * ((AgSigma ** 12) / (R ** 13) - 0.5 * (AgSigma ** 6) / (R ** 7)))

        VelX = -(Positions[0][0] - Positions[i][0]) * Force / R
        VelY = -(Positions[0][1] - Positions[i][1]) * Force / R
        VelZ = -(Positions[0][2] - Positions[i][2]) * Force / R

        forces[0][0] += VelX
        forces[0][1] += VelY
        forces[0][2] += VelZ

    '''Вычисление силы из потенциала Морзе'''
    for i in range(1, numOfAtoms):
        for j in range(i+1,numOfAtoms):

            R = math.sqrt((Positions[i][0] - Positions[j][0]) ** 2 +
                          (Positions[i][1] - Positions[j][1]) ** 2 +
                          (Positions[i][2] - Positions[j][2]) ** 2)
            Force = AlEps*math.exp(-2*AlAlpha*R)*(AlAlpha*math.exp(AlAlpha*R)-4*AlAlpha)


            VelX = -(Positions[i][0] - Positions[j][0]) * Force / R
            VelY = -(Positions[i][1] - Positions[j][1]) * Force / R
            VelZ = -(Positions[i][2] - Positions[j][2]) * Force / R

            forces[i][0] += VelX
            forces[i][1] += VelY
            forces[i][2] += VelZ

            forces[j][0] -= VelX
            forces[j][1] -= VelY
            forces[j][2] -= VelZ

    return forces

'''Применение сил к молекулам по второму закону Ньютона'''
def applyForces(Positions, Velocities, forces, Masses, TimeStep):
    print("Vel1: ", Velocities* TimeStep)
    Positions += Velocities * TimeStep
    Velocities += forces * TimeStep / Masses[np.newaxis].T
    print("Vel2: ", Velocities* TimeStep)

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
    AgVelocity = np.array([0.0, 0.0, 0.0])

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

        forces = calcForces(Velocities, Positions) # Вычисление сил
        applyForces(Positions, Velocities, forces, Masses, TimeStep) # Интегрирование позиций и скоростей
        print("Pos:", Positions)
        Dump.writeOutput(OutputFileName, AlNumOfAtoms + 1, step, Borders,
                         radius=Radius, pos=Positions, v=Velocities)  # Вывод данных в .dump файл
