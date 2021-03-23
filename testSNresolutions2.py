import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

if __name__ == '__main__':
    radii = [r for r in range(4,41)]
    #with open('resolutionOscillationsForTinyCircles2.csv','w',newline='') as file:
    with open('resolutionTestExpanseMar21.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['readshock2/rs2p5_1.csv'])
        writer.writerow(['inclusive dist. mar 21. trying to just see how runtimes on expanse look'])
        writer.writerow(['rad','S','N','S/N','runtime'])

        for radius in radii:
            coll = exclVol.PolycrystalGrid('crystals/readshock2/rs2p5_1.csv',rad=radius,usePsi6=False)
            print(coll.beadRad)
            (S,numParts,Si,time) = coll.entropyParallel(40)
            writer.writerow([radius,S,numParts,S/numParts,time])
