import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)



coll = exclVol.PolycrystalGrid('tinyCirc1_10.csv',rad=8,usePsi6=True)

numSnow = 0
for i in range(len(coll.particleCenters)):
    if coll.psi6s[i] >= 0.98:
        numSnow += 1
print(numSnow)
print(len(coll.particleCenters))

(S,Sbead,numParts,time) = coll.entropyHisto()
print("time: "+str(time))
print("S/N: "+str(S/numParts))
print("S_snowflake: "+str(np.log( len(coll.snowflakeShape)/len(coll.beadShape) )))

Smain = []
Ssnow = []
for (Si, psi6) in Sbead:
    if psi6 >= 0.98:
        Ssnow.append(Si)
    else:
        Smain.append(Si)

print(len(Ssnow))
print(len(Smain))

plt.hist([Smain,Ssnow],stacked=True)
plt.show()