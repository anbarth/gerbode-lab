import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)



coll = exclVol.PolycrystalGrid('tinyTest3_10.csv',rad=10,usePsi6=True)
coll.showGrid()