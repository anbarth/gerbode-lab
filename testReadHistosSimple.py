import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)


cutoff = 0.98
dir1 = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\rs2_Si_5px'
dir2 = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\rs2_Si_35px'


for filename in os.listdir(dir1):
    print(filename)

    Sbead = []
    with open('rs2_Si_5px/'+filename) as file:
        reader = csv.reader(file)
        #next(reader) # these files do not have Ssnow at the top
        for row in reader:
            Sbead.append(float(row[0]))

    plt.clf()
    plt.xlim(-4.5,-1)
    plt.ylim(0,600)
    plt.hist(Sbead)
    plt.vlines(np.mean(Sbead),0,25)
    plt.savefig('rs2_5px_35px_compare_bounded/'+filename[0:-4]+'_5px_histo.png')

for filename in os.listdir(dir2):
    print(filename)

    Sbead = []
    with open('rs2_Si_35px/'+filename) as file:
        reader = csv.reader(file)
        next(reader) # these files have Ssnow at the top
        for row in reader:
            Sbead.append(float(row[0]))
    
    plt.clf()
    plt.xlim(-4.5,-1)
    plt.ylim(0,600)
    plt.hist(Sbead)
    plt.vlines(np.mean(Sbead),0,50)
    plt.savefig('rs2_5px_35px_compare_bounded/'+filename[0:-4]+'_35px_histo.png')