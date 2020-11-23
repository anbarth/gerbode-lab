import csv
import numpy as np

fname = 'small.csv'

beadRad = 4
spacing = 2

with open(fname,'w',newline='') as csvFile:
    writer = csv.writer(csvFile,delimiter=',')
    x = beadRad
    for i in range(4):
        if i % 2 == 0:
            y = beadRad
        else:
            y = beadRad+(beadRad*2 + spacing)/2
        for j in range(4):
            writer.writerow([x,y])
            y += beadRad*2 + spacing
            #y += (beadRad*2 + spacing)/2
        x += (beadRad*2 + spacing)*np.sqrt(3)/2
        