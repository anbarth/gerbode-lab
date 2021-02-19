import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

coll = exclVol.PolycrystalGrid('annasCoolTest10.csv',resolution=9/4.5,usePsi6=True)
print(coll.entropy())

coll2 = exclVol.PolycrystalGrid('annasCoolTest10.csv',resolution=9/4.5,usePsi6=False)
print(coll2.entropy())
