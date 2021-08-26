import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)


def makeEnergyBarrier(simName):
    dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName

    with open("../grainsplits/"+simName+"/energyCurve.csv",'w',newline='') as dataOut:
        writer = csv.writer(dataOut)
        writer.writerow(['file','timestep','S','N','F'])
        step = 0
        for filename in os.listdir(dir):  
            print(filename)
            coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride=True)
            S = coll.entropy(imgFile="../grainsplits/"+simName+"/snowflakes/"+str(filename[0:-4])+".png")
            N = coll.numParts()
            F = -1*S
            writer.writerow([filename,step,S,N,F])
            step = step+1


#[[10, 3], [10, 5], [10, 7], 
#      [15, 3], [15, 5], [15, 7], [15, 9], [15, 12], 
#      [20, 3], [20, 5], [20, 7], [20, 9], [20, 12], [20, 15], 
#      [25, 3], [25, 5], [25, 7], [25, 9], [25, 12], [25, 15], 
Dds =  [[30, 3], [30, 5], [30, 7], [30, 9], [30, 12], [30, 15]]

for (D,d) in Dds:
    simName = "aug2/D"+str(D)+"/d"+str(d)
    makeEnergyBarrier(simName)