import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


if __name__ == '__main__':
    for filename in ['tinyCircle1_8.csv','tinyCircle1_10.csv']:
        coll = exclVol.PolycrystalGrid(filename,rad=35,usePsi6=True)
        (S,Sbead,numParts,time) = coll.entropyParallelHisto(40)
        print(len(coll.particleCenters))
        print("time: "+str(time))
        with open(filename[0:-4]+'_Sbead.csv','w',newline='') as file:
            writer = csv.writer(file)
            # first line: S_snowflake
            writer.writerow([np.log( len(coll.snowflakeShape)/len(coll.beadShape) )])
            # remaining lines: (S_i with NO shortcut,|psi6|)
            for (Si,psi6) in Sbead:
                writer.writerow([Si,psi6])