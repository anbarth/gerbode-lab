import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\readshock'
with open('edgesOrNoEdges.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['file','S_edges','N_edges','S/N_edges','S_noedges','N_noedges','S/N_noedges'])
    for filename in os.listdir(dir):
        if (filename.startswith("readshock1") or filename.startswith("readshock2")) and filename.endswith(".csv"):
            print(filename)
            coll = PolycrystalGrid('readshock/'+filename,resolution=9/5)
            [S1,N1,runtime] = coll.entropy()
            [S2,N2] = coll.entropyNoEdgeBeads()
            writer.writerow([filename,S1,N1,S1/N1,S2,N2,S2/N2])