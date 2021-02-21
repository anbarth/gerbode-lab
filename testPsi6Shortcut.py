import csv
import exclVol
import importlib
import os
importlib.reload(exclVol)

if __name__ == '__main__':

    dir = r'C:\Users\GerbodeLab\Documents\banana\gerbode-lab\readshock2'
    with open('psi6shortcutTest2.csv','w',newline='') as dataOut:
        writer = csv.writer(dataOut)
        writer.writerow(['file','S_nocut','N_nocut','S/N_nocut','t_nocut','S_cut','N_cut','S/N_cut','t_cut'])
        #for filename in os.listdir(dir):
        for filename in ['rs10_1.csv']:
            print(filename)
            coll = exclVol.PolycrystalGrid('readshock2/'+filename,resolution=35/5)
            [S1,N1,time1] = coll.entropyParallel(40)
            coll2 = exclVol.PolycrystalGrid('readshock2/'+filename,resolution=35/5,usePsi6=True)
            [S2,N2,time2] = coll2.entropyParallel(40)
            writer.writerow([filename,S1,N1,S1/N1,time1,S2,N2,S2/N2,time2])