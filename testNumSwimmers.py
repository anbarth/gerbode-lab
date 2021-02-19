import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)


if __name__ == '__main__':
    resolutions = [4/5,1,6/5,7/5,8/5,9/5,2,11/5,12/5,13/5,14/5,3,16/5,17/5,18/5,19/5,20/5]
    #resolutions = [5/10,6/10,7/10,8/10,9/10,1,11/10,12/10,13/10]
    #resolutions = [3/4,1,5/4,6/4,7/4,8/4,9/4,10/4,11/4,12/4,13/4,14/4,15/4,16/4,17/4]
    with open('megaResolutionData.csv','w',newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['beadRad','S_80','t_80'])
        for res in resolutions:
            coll = PolycrystalGrid('readshock/readshock4_7p5.csv',resolution=res)
            print(coll.beadRad)
            [S, N, t] = coll.entropyParallel(80)
            #[S10, N10, t10] = coll.entropyParallel(10)
            #[S20, N20, t20] = coll.entropyParallel(20)
            #[S40, N40, t40] = coll.entropyParallel(40)
            writer.writerow([coll.beadRad,S,t])