import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)


if __name__ == '__main__':
    dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\crystals\readshock2'
    with open('readshock2_poly_45px_mar26.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['readshock2 data w freeSpacePoly, inclusive dist, 40 worker processes. what will the trend be??'])
        writer.writerow(['file','S','N','S/N','time'])
        for filename in os.listdir(dir):
            print("==================")
            print(filename)
            coll = exclVol.PolycrystalGrid('crystals/readshock2/'+filename,rad=45,useNeighbs=True)
            (S,numParts,Sbead,time) = coll.entropyParallel(40,makeImg=True,poly=True)
            print(S,numParts,time)
            writer.writerow([filename,S,numParts,S/numParts,time])
            # expected runtime 9 hours