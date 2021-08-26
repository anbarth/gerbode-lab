import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

R_grain = 60
tightness = 1.05
rad = 5

xmax = 2*R_grain+200
ymax = 2*R_grain+200
buffer = rad*tightness

theta = np.linspace(0,2*np.pi,25)
x = xmax/2 + (R_grain+buffer)*cos(theta)
y = ymax/2 + (R_grain+buffer)*sin(theta)

angleNames = ['5','6p25','7p5','8p75','10','15','20','25']
simNames = ['r60/theta'+aName for aName in angleNames]

for simName in simNames:
    "crystals/"+simName+"/0000.csv"
    windowPath = path.Path(self.windowVertices)
    inWindow = windowPath.contains_points(self.particleCenters)