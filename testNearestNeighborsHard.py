import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
importlib.reload(exclVol)

coll = PolycrystalGrid('nearestNeighborsHard.csv')
space = sorted(coll.freeSpace((20,15)))
with open('nearestNeighborConfigs2.csv','w',newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['num beads','num configs','runtime'])
    for n in range(1,6):
        print(n)
        tic = time.time()
        configs = coll.numConfigs(n,space)
        toc = time.time()
        writer.writerow([n,configs,toc-tic])