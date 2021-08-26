import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib import path
import polycrystal
import importlib
import os
importlib.reload(polycrystal)


#dir = "C:\\Users\\anna2\\OneDrive\\Documents\\Gerbode\\python\\crystals\\flatGBs2"

#newRadii = [0,4.9935,4.9995,4.994,4.997,4.9967,5,4.9947,4.9977,4.9873,4.992,4.996] # flatGBs2
#newRadii = [0,4.9913,5,4.9896,4.9947,4.9953,4.9983,4.9939,4.9928,4.983,4.9941,4.995] # flatGBs3
newRadii = [4.9841,4.986,5,4.9886,4.9869,4.9942,4.9843,4.9843,4.9941,4.9894,4.9936,4.9918] # flatGBs4


with open("../flatGBs4/flatGBEnergies_blah.csv",'w',newline='') as dataOut:
    writer = csv.writer(dataOut)
    writer.writerow(['file','miso','N','S','F','S_corrected','F_corrected'])
    step = 0
    for i in range(1,2):
        miso = 2.5*i
        misoStr = str(miso)
        if i % 2 == 0:
            misoStr = str(int(miso))
        filename = misoStr+".csv"
        print(filename)

        coll = polycrystal.Polycrystal("crystals/flatGBs4/"+filename,windowOverride=False,radius=5)
        S = coll.entropy(imgFile="../flatGBs4/snowflakes/"+str(filename[0:-4])+"_snowflakes.png")

        myRad = newRadii[i-1]
        coll = polycrystal.Polycrystal("crystals/flatGBs4/"+filename,windowOverride=False,radius=myRad)
        S_corr = coll.entropy(imgFile="../flatGBs4/snowflakes_new/"+str(filename[0:-4])+"_snowflakes.png")
        N = coll.numParts()
        writer.writerow([filename,miso,N,S,-S,S_corr,-S_corr])