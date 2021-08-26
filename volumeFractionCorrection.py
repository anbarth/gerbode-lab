import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

simName = 'r50theta20_lloyd/therm01'
dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName

radInit = 5
radFinal = 4.99
radStep = (radFinal-radInit)/10

with open("deleteme.csv",'w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['file','timestep','bead rad','S','N','F'])
    step = 0
    S_0 = 0
    for filename in os.listdir(dir):  
        print(filename)
        myRad = radInit + step*radStep
        coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride=True,radius=myRad)
        [S,N,t] = coll.entropy(imgFile=str(filename[0:-4])+".png")
        if step == 0:
            S_0 = S
        F = -1*S
        writer.writerow([filename,step,myRad,S,N,F])
        step = step+1