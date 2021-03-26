import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


if __name__ == '__main__':
    radii = [40,45]
    with open('polyScalingMar25.csv','w',newline='') as file:
        writer = csv.writer(file)
        crysfilename = 'crystals/readshock2/rs7p5_1.csv'
        writer.writerow([crysfilename])
        writer.writerow(['can we take resolution even higher?? inclusive distance. 40 worker processes.'])
        writer.writerow(['rad','N','S_old','runtime','S_poly','runtime'])
        
        for radius in radii:
            coll = exclVol.PolycrystalGrid(crysfilename,rad=radius,useNeighbs=True)
            print(coll.beadRad)
            (Sold,N1,Sbead1,told) =   coll.entropyParallel(40,makeImg=True,poly=False)
            print("done with old freeSpace")
            (Spoly,N2,Sbead2,tpoly) = coll.entropyParallel(40,makeImg=True,poly=True)
            print(radius,N1,Sold,told,Spoly,tpoly)
            writer.writerow([radius,N1,Sold,told,Spoly,tpoly])