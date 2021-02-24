import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


if __name__ == '__main__':
    coll = exclVol.PolycrystalGrid('tinyCircle1_8.csv',resolution=5,usePsi6=True)
    coll.showGrid()
