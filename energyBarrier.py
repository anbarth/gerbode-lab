import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

simName = 'jul4/r75/theta20/sim07'
dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName

with open('../grainsplits/'+simName+'/barrier.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    #writer.writerow(['peanut exp entropy values found 6/16/2021, frame 26-frame 27 with intermediate timesteps found via moveparts.m'])
    writer.writerow(['file','timestep','S','N','F'])
    step = 0
    S_0 = 0
    for filename in os.listdir(dir):  
        print(filename)
        coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride=True)
        #coll.show()
        [S,N,t] = coll.entropy(imgFile="../grainsplits/"+simName+"/snowflakes/"+str(filename[0:-4])+".png")
        if step == 0:
            S_0 = S
        F = -1*S + S_0
        writer.writerow([filename,step,S,N,F])
        step = step+1