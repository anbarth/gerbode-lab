import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\crystals\peanut'
radii = [10,10.25,10.5,10.75,11]
radStr = ["01_10","02_1025","03_105","04_1075","05_11"]
with open('deleteme.csv','w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    #writer.writerow(['peanut exp entropy values found 6/16/2021, frame 26-frame 27 with intermediate timesteps found via moveparts.m'])
    writer.writerow(['radius','file','S','N'])
    for i in range(len(radii)):
        rad = radii[i]
        myRadStr = radStr[i]
        for filename in os.listdir(dir):  
            print(rad,filename)      
            coll = polycrystal.Polycrystal("crystals/peanut/"+filename,radius=rad)#,windowOverride=True)
            [S,N,t] = coll.entropyDraw(imgFile="../peanut_experiment/"+myRadStr+"_"+str(filename[0:-4])+".png")
            writer.writerow([coll.beadRad,filename,S,N])