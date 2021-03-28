import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


if __name__ == '__main__':
    radii = [r for r in range(5,48)]
    with open('polyScalingMar25.csv','w',newline='') as file:
        writer = csv.writer(file)
        crysfilename = 'crystals/readshock2/rs2p5_1.csv'
        writer.writerow([crysfilename])
        writer.writerow(['freeSpacePoly, inclusive distance. how does S/N plateau? 40 worker processes. DrSpeedy'])
        writer.writerow(['rad','S_poly','N','runtime'])
        
        for radius in radii:
            coll = exclVol.PolycrystalGrid(crysfilename,rad=radius,useNeighbs=True)
            print(coll.beadRad)
            (Spoly,N,Sbead,t) = coll.entropyParallel(40,makeImg=False,poly=True)
            print(radius,Spoly,N,t)
            writer.writerow([radius,Spoly,N,t])