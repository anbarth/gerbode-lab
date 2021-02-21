import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)



for fname in ['annasCoolTest10.csv']:
    coll = exclVol.PolycrystalGrid(fname,resolution=5/4.5,usePsi6=True)
    (S,Smain,Ssnow,numParts,time) = coll.entropyHisto()
    with open(fname[0:-4]+'_Smain.csv','w',newline='') as file:
        writer = csv.writer(file)
        for ent in Smain:
            writer.writerow([ent])
    with open(fname[0:-4]+'_Ssnow.csv','w',newline='') as file:
        writer = csv.writer(file)
        for ent in Ssnow:
            writer.writerow([ent])
    #plt.hist([Smain,Ssnow],stacked=True)
    #plt.show()

#plt.hist(coll.psi6s,range=[0.9,1])
#plt.show()