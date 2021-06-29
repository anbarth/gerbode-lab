import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

simNames = ['r50_theta'+str(5*x) for x in range(1,6)]
with open('../grainsplits/rs50barriers.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['sim name','file','timestep','S','N'])

    for simName in simNames:
        print("=========================")
        print(simName)
        dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName
        step = 0
        for filename in os.listdir(dir):  
            print(filename)      
            coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride=True)
            #coll.show()
            [S,N,t] = coll.entropy(imgFile="../grainsplits/"+simName+"/snowflakes/"+str(filename[0:-4])+".png")
            writer.writerow([simName,filename,step,S,N])
            step = step+1