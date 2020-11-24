import csv
import numpy as np

fname = 'nearestNeighbors.csv'

beadRad = 4
spacing = 2
# TODO you gotta manually go into the csv file and add a row at the top for [grid X, grid Y]

with open(fname,'w',newline='') as csvFile:
    writer = csv.writer(csvFile,delimiter=',')
    x = beadRad
    for i in range(5):
        if i % 2 == 0:
            y = beadRad
        else:
            y = beadRad+(beadRad*2 + spacing)/2
        for j in range(5):
            writer.writerow([x,y])
            y += beadRad*2 + spacing
            #y += (beadRad*2 + spacing)/2
        x += (beadRad*2 + spacing)*np.sqrt(3)/2
        