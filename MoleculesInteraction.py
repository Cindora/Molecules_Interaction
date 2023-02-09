import Dump
import numpy as np
import os


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

    # print(AlPositions)

    for i in range(dimension):  # Подгонка значений относительно границ области
        AlPositions[:, i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * AlPositions[:, i]
        AgPosition[i] = Borders[i][0] + (Borders[i][1] - Borders[i][0]) * AgPosition[i]

    # print("Pos: ", AlPositions, "\n", AgPosition)
    # print("Vel: ", AlVelocity, "\n", AgVelocity)
    Positions = np.append([AgPosition], AlPositions, axis=0)
    Velocities = np.append([AgVelocity], AlVelocity, axis=0)
    Radius = np.append([AgRadius], AlRadius)

    # print(Positions, "\n", Velocities, "\n" ,Radius)

    if os.path.exists(OutputFileName):
        os.remove(OutputFileName)

    step = 0

    Dump.writeOutput(OutputFileName, AlNumOfAtoms + 1, step, Borders,
                     radius=Radius, pos=Positions, v=Velocities)  # Вывод данных в .dump файл
