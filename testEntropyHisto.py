import numpy as np
import csv
import matplotlib.pyplot as plt
import polycrystal
import importlib
import os
importlib.reload(polycrystal)

if __name__ == '__main__':
    dir = r'C:\Users\anna2\OneDrive\Documents\Gerbode\python\crystals\readshock2'
    #dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock'
    with open('deleteme.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['readshock2 data, finally grid-free! see also associated images'])
        #writer.writerow(['readshock2 data at 5px, see also associated images!'])
        writer.writerow(['file','S','N','S/N','time'])
        for filename in os.listdir(dir):
            print("==================")
            print(filename)
            coll = polycrystal.Polycrystal('crystals/readshock2/'+filename)
            (S,numParts,time) = coll.entropyDraw()
            #print("time: "+str(time))
            #writer.writerow([filename,S,numParts,S/numParts,time])
    