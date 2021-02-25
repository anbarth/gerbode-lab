import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)


cutoff = 0.98
#dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\Sbead_data'
with open('whySnowflakeBadTestNew.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow([cutoff])
    writer.writerow(['file','Ntot','Nordered','avg(Sdiso)','avg(Sordered)','Ssnow','stddev(Sdiso)','stddev(Sord)','S/N_actual','S/N_cutoff'])
    #for filename in os.listdir(dir):
    for filename in ['tinyCircle1_8_Sbead.csv','tinyCircle1_10_Sbead.csv','tinyCircle2_8_Sbead.csv','tinyCircle2_10_Sbead.csv']:
        print(filename)

        Sbead = []
        Sdiso = []
        Sord = []
        #with open('Sbead_data/'+filename) as file:
        with open(filename) as file:
            reader = csv.reader(file)
            Ssnow = float(next(reader)[0])
            for row in reader:
                Sbead.append([float(row[0]),float(row[1])])
        
        for (S,psi6) in Sbead:
            if psi6 >= cutoff:
                Sord.append(S)
            else:
                Sdiso.append(S)
        
        Sactual = np.mean(Sdiso+Sord)
        Scutoff = (len(Sdiso)*np.mean(Sdiso) + len(Sord)*Ssnow)/len(Sbead)


        writer.writerow([filename,len(Sbead),len(Sord), \
                        np.mean(Sdiso),np.mean(Sord), Ssnow, \
                        np.std(Sdiso),np.std(Sord), Sactual, Scutoff])
        
        plt.clf()
        #plt.hist(Sord)
        plt.hist(np.array([Sdiso,Sord],dtype=object),stacked=True)
        #plt.vlines(Ssnow,0,50)
        plt.savefig(filename[0:-4]+'_histo.png')