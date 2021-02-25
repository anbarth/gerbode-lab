import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)

if __name__ == '__main__':
    #dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock2'
    #for filename in os.listdir(dir):
    for filename in ['tinyCircle1_8.csv','tinyCircle1_10.csv','tinyCircle2_8.csv','tinyCircle2_10.csv']:
        print("==================")
        print(filename)
        #coll = exclVol.PolycrystalGrid('readshock2/'+filename,resolution=35/5,usePsi6=True)
        coll = exclVol.PolycrystalGrid(filename,rad=35,usePsi6=True)
        print("n_bead: "+str(len(coll.beadShape)))
        print("n_snow: "+str(len(coll.snowflakeShape)))
        (S,Sbead,numParts,time) = coll.entropyParallelHisto(40)
        print("time: "+str(time))
        with open(filename[0:-4]+'_Sbead.csv','w',newline='') as file:
            writer = csv.writer(file)
            # first line: S_snowflake
            writer.writerow([np.log( len(coll.snowflakeShape)/len(coll.beadShape) )])
            # remaining lines: (S_i with NO shortcut,|psi6|)
            for (Si,psi6) in Sbead:
                writer.writerow([Si,psi6])