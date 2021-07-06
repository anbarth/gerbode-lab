import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

#simNames = ['r50_theta'+str(5*x) for x in range(1,6)]
angleNames = ['5','6p25','7p5','8p75','10','15','20','25']
#angleNames = ['8p75','10','15','20','25']
simNames = ['r60/theta'+aName for aName in angleNames]

with open('../grainsplits/r60barriersWithCircularWindow.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['sim name','file','timestep','S','N','F'])

    for simName in simNames:
        print("=========================")
        print(simName)
        dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\"+simName
        step = 0
        S_0 = 0
        os.mkdir("../grainsplits/"+simName+"/snowflakesCircularWindowConstParticles")
        for filename in os.listdir(dir):  
            print(filename)  
            coll = polycrystal.Polycrystal("crystals/"+simName+"/"+filename,windowOverride="../grainsplits/"+simName+"/partsInCircle.csv")
            #coll.show()
            [S,N,t] = coll.entropy(imgFile="../grainsplits/"+simName+"/snowflakesCircularWindowConstParticles/"+str(filename[0:-4])+".png")
            if step == 0:
                S_0 = S
            F = -1*S + S_0
            writer.writerow([simName,filename,step,S,N,F])
            step = step+1