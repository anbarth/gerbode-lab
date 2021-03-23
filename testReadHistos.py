import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
import PIL
importlib.reload(exclVol)


dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\Sbeads\rs2_35px_excl_mar17'
with open('deleteme.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    #writer.writerow(['rs2 entropy values found mar 10 and 11. not necessarily the same as the last time i did it bc remember that little bug i fixed with how freespace starts out'])
    #writer.writerow(['file','S','N','S/N'])
    for filename in os.listdir(dir):
        Sbead = []
        with open('Sbeads/rs2_35px_excl_mar17/'+filename) as file:
            reader = csv.reader(file)
            for row in reader:
                Sbead.append(float(row[0]))

        S = np.mean(Sbead)
        N = len(Sbead)

        #writer.writerow([filename,S,N,S/N])

        plt.clf()
        plt.hist(Sbead,bins=40)
        plt.vlines(S,0,50)
        plt.xlim(-4.5,-1)
        plt.ylim(0,170)
        plt.savefig(filename[0:-4]+'_35px_EXCL_histo.png')