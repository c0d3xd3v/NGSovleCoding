import math
import ngsolve

from pde.helmholtz import *


def solveAcousticBoundaryValue(air_fes, n, g,roh, boundaries):
    gfu = GridFunction(air_fes, name='gfu')
    gfu.Set(roh*g*n, definedon=boundaries)
    gfu = solveHelmholtz(air_fes, gfu, roh)
    return gfu

def rohForAir():
    # The air is the acoustic medium at the room temperature with a 
    # density of 1.21 kgm−3 and sound speed of 340 ms−1 at a 
    # reference pressure of 20 μPa
    # https://www.mdpi.com/2624-599X/5/1/10
    density_air = 1.013 # bar
    compressibility_air = 101.0 # bar
    f = 540.0
    roh = 2.0 * math.pi * f * sqrt(compressibility_air * density_air)
    roh = roh*0.0000015*0.125
    return roh
