import numpy as np
import random
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\Sbead_data'
for filename in os.listdir(dir):
#for filename in ['rs10_1_Sbead.csv']:
    Sbead = []
    Smain = []
    Spsihigh = []
    with open('Sbead_data/'+filename) as file:
        reader = csv.reader(file)
        Ssnow = float(next(reader)[0])
        for row in reader:
            Sbead.append([float(row[0]),float(row[1])])
    
    for (S,psi6) in Sbead:
        if psi6 >= 0.98:
            Spsihigh.append(S)
        else:
            Smain.append(S)
    
    Ssnowrunning = []
    j = random.randrange(0,len(Spsihigh))
    Ssnowrunning.append(Spsihigh[j])
    for i in range(399):
        j = random.randrange(0,len(Spsihigh))
        Ssnowrunning.append( (Ssnowrunning[i]*(i+1)+Spsihigh[j])/(i+2)  )
    #print(len(Spsihigh))
    #for i in range(1,len(Spsihigh)+1):
    #    Ssnowrunning.append(np.mean(Spsihigh[0:i]))


    plt.clf()
    #plt.xlim(0,200)
    plt.hlines(np.mean(Spsihigh),1,400)
    plt.hlines(np.mean(Spsihigh)*1.02,1,400,ls='--')
    plt.hlines(np.mean(Spsihigh)*0.98,1,400,ls='--')
    #plt.plot([i for i in range(1,201)],Ssnowrunning)
    plt.plot([i for i in range(len(Ssnowrunning))],Ssnowrunning)
    plt.savefig(filename[0:-4]+'_runningSnowRando.png')

