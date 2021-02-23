import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

Smain = []
with open('rs10_1_Smain.csv') as file:
    reader = csv.reader(file)
    for row in reader:
        Smain.append(float(row[0]))

Ssnow = []
with open('rs10_1_Ssnow.csv') as file:
    reader = csv.reader(file)
    for row in reader:
        Ssnow.append(float(row[0]))

#plt.hist([Smain,Ssnow],stacked=True,range=[-4,-2])
#plt.hist(Ssnow,stacked=True,range=[-4,-2])
plt.show()