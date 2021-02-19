import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


resolutions = [15,16,17,18,19,20]
with open('areaFractions4.csv','w',newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['beadRad','occupied px','px in square','fraction'])
    for res in resolutions:
        coll = PolycrystalGrid('oneBead.csv',resolution=res)
        occ = len(coll.occupiedPx)
        square = 4*coll.beadRad*coll.beadRad
        writer.writerow([coll.beadRad,occ,square,occ/square])