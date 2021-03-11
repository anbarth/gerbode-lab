import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

if __name__ == '__main__':
    radii = [r for r in range(4,31)]+[32,35,37,40]
    #with open('resolutionOscillationsForTinyCircles2.csv','w',newline='') as file:
    with open('megaResolutionTestExclusive.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow('readshock/readshock3_15.csv')
        writer.writerow(['rad','S','N','S/N','runtime'])

        for radius in radii:
            coll = exclVol.PolycrystalGrid('readshock/readshock3_15.csv',rad=radius,usePsi6=False)
            print(coll.beadRad)
            (S,numParts,time,Si) = coll.entropyParallel(40)
            writer.writerow([radius,S,numParts,S/numParts,time])