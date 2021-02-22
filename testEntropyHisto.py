import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)

if __name__ == '__main__':
    dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock2'
    #for filename in os.listdir(dir):
    for filename in ['rs10_1.csv']:
        print(filename)
        coll = exclVol.PolycrystalGrid('readshock2/'+filename,resolution=5/5,usePsi6=True)
        (S,Smain,Ssnow,numParts,time) = coll.entropyHisto()
        print(time)
        '''with open(filename[0:-4]+'_Smain.csv','w',newline='') as file:
            writer = csv.writer(file)
            for ent in Smain:
                writer.writerow([ent])
        with open(filename[0:-4]+'_Ssnow.csv','w',newline='') as file:
            writer = csv.writer(file)
            for ent in Ssnow:
                writer.writerow([ent])'''
        plt.hist([Smain,Ssnow],stacked=True)
        plt.show()
        plt.clf()

        (S,Smain,Ssnow,numParts,time) = coll.entropyParallelHisto(40)
        print(time)
        plt.hist([Smain,Ssnow],stacked=True)
        plt.show()
        plt.clf()

    #plt.hist(coll.psi6s,range=[0.9,1])
    #plt.show()