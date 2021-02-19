import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

if __name__ == '__main__':
    coll = exclVol.PolycrystalGrid('readshock2/rs5_1.csv')
    print(coll.entropy())
    print(coll.entropyParallel(40))

    coll2 = exclVol.PolycrystalGrid('readshock2/rs5_1.csv',usePsi6=True)
    print(coll2.entropy())
    print(coll2.entropyParallel(40))

    #dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock2'
    '''with open('psi6shortcutTest.csv','w',newline='') as dataOut:
        writer = csv.writer(dataOut)
        writer.writerow(['file','S_hard','N_hard','S/N_edges','S_noedges','N_noedges','S/N_noedges'])
        for filename in os.listdir(dir):
            if (filename.startswith("readshock1") or filename.startswith("readshock2")) and filename.endswith(".csv"):
                print(filename)
                coll = PolycrystalGrid('readshock/'+filename,resolution=9/5)
                [S1,N1,runtime] = coll.entropy()
                [S2,N2] = coll.entropyNoEdgeBeads()
                writer.writerow([filename,S1,N1,S1/N1,S2,N2,S2/N2])'''