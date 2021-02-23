import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)



coll = exclVol.PolycrystalGrid('readshock2/'+'rs10_1.csv',resolution=35/5,usePsi6=True)
#coll.showGrid()