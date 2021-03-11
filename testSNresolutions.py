import numpy as np
import csv
import matplotlib.pyplot as plt
import exclVol
import importlib
importlib.reload(exclVol)

if __name__ == '__main__':
    radii = [r for r in range(4,39)]
    with open('resolutionOscillationsForTinyCircles2.csv','w',newline='') as file:
    #with open('readshock/readshock3_15.csv') as file:
        writer = csv.writer(file)
        writer.writerow(['filename','rad','S','N','S/N','runtime'])
        for crysfilename in ['tinyCircle1_8.csv','tinyCircle1_10.csv','tinyCircle2_8.csv','tinyCircle2_10.csv']:
            print(crysfilename)

            for radius in radii:
                coll = exclVol.PolycrystalGrid(crysfilename,rad=radius,usePsi6=True)
                print(coll.beadRad)
                (S,Sbead,numParts,time) = coll.entropyParallelHisto(40)

                with open(crysfilename[0:-4]+'_'+str(radius)+'_Sbead.csv','w',newline='') as sbeadfile:
                    writer2 = csv.writer(sbeadfile)
                    # first line: S_snowflake
                    writer2.writerow([np.log( len(coll.snowflakeShape)/len(coll.beadShape) )])
                    # remaining lines: (S_i with NO shortcut,|psi6|)
                    for (Si,psi6) in Sbead:
                        writer2.writerow([Si,psi6])

                #[S2, N2, runtime2] = coll.entropyParallel()
                writer.writerow([crysfilename,radius,S,numParts,S/numParts,time])