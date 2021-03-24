import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

if __name__ == '__main__':
    radii = [5,10,15,20,25,30,35]
    with open('freeSpaceOldVsPolyMar24.csv','w',newline='') as file:
        writer = csv.writer(file)
        crysfilename = 'crystals/tinyCircle1_8.csv'
        writer.writerow([crysfilename])
        writer.writerow(['supplement for freeSpaceOldVsPolyMar23, now with exclusive dist'])
        writer.writerow(['rad','N','S_old','runtime'])
        
        for radius in radii:
            coll = exclVol.PolycrystalGrid(crysfilename,rad=radius,useNeighbs=True)
            print(coll.beadRad)
            (Sold,N1,Sbead1,told) =   coll.entropyParallel(40,makeImg=True,poly=False)
            #(Spoly,N2,Sbead2,tpoly) = coll.entropyParallel(40,makeImg=True,poly=True)
            #print(N1==N2) # true probly
            print(radius,N1,Sold,told)
            writer.writerow([radius,N1,Sold,told])