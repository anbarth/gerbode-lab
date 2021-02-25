import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
import os
importlib.reload(exclVol)



dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\Sbead_data'
with open('whySnowflakeBad99.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['file','Ntot','Npsihigh','avg(Smain)','avg(Spsihigh)','Ssnow','stddev(Smain)','stddev(Spsihigh)'])
    for filename in os.listdir(dir):
        Sbead = []
        Smain = []
        Spsihigh = []
        with open('Sbead_data/'+filename) as file:
            reader = csv.reader(file)
            Ssnow = float(next(reader)[0])
            for row in reader:
                Sbead.append([float(row[0]),float(row[1])])
        print(len(Sbead))
        
        for (S,psi6) in Sbead:
            if psi6 >= 0.99:
                Spsihigh.append(S)
            else:
                Smain.append(S)
        
        writer.writerow([filename,len(Sbead),len(Spsihigh), \
                        np.mean(Smain),np.mean(Spsihigh), Ssnow, \
                        np.std(Smain),np.std(Spsihigh)])
        
        plt.clf()
        plt.hist(Spsihigh)
        #plt.hist(np.array([Smain,Spsihigh],dtype=object),stacked=True)
        #plt.vlines(Ssnow,0,50)
        plt.savefig(filename[0:-4]+'_histoSnow99.png')


        


#plt.hist(Spsihigh)#,range=[-4,-2])
#plt.show()