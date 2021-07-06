import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)



simNames = ['jul4/r75/theta5/sim'+f'{i:02}' for i in range(1,11)]


with open('delteme.csv','w',newline='') as dataOut:
#with open('../grainsplits/jul4/r75theta20_10sims.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    #writer.writerow(['sim','file','timestep','S','N','F'])

    for simName in simNames:
        print("=========================")
        print(simName)
        dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName
        step = 0
        S_0 = 0
        #os.mkdir("../grainsplits/jul4/"+simName+"/snowflakesCircularWindowConstParticles")
        for filename in os.listdir(dir):  
            print(filename)  
            coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride=True)
            #coll.show()
            [S,N,t] = coll.entropy(imgFile="../grainsplits/"+simName+"/snowflakes/"+str(filename[0:-4])+".png")
            if step == 0:
                S_0 = S
            F = -1*S + S_0
            #writer.writerow([simName,filename,step,S,N,F])
            step = step+1