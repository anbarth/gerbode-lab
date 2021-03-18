import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)

if __name__ == '__main__':
    dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\crystals\readshock2'
    #dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock'
    with open('readshock2_35px_excl_Entropiesmarch17.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['readshock2 data at 35px, see also associated images! this time EXCLUSIVE dist!'])
        #writer.writerow(['readshock2 data at 5px, see also associated images!'])
        writer.writerow(['file','S','N','S/N','time'])
        for filename in os.listdir(dir):
            print("==================")
            print(filename)
            coll = exclVol.PolycrystalGrid('crystals/readshock2/'+filename,rad=35,usePsi6=True)
            #coll = exclVol.PolycrystalGrid(filename,rad=35)
            print("n_bead: "+str(len(coll.beadShape)))
            (S,numParts,Sbead,time) = coll.entropyParallel(40,makeImg=True)
            print("time: "+str(time))
            writer.writerow([filename,S,numParts,S/numParts])
    