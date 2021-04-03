import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import bisect
import random
import operator
import os
import multiprocessing as mp
import PIL
import myGeo
import importlib
from matplotlib import path
importlib.reload(myGeo)

class Polycrystal:
    # class fields:
    #   xmin: minimum x value
    #   ymin: minimum y value TODO explain these better lol
    #   windowSize: a list that goes [width,height]. note that 0 is included and width and height are excluded.
    #   beadRad: bead radius in px
    #   particleCenters: 

    # constructor takes the input csv file
    def __init__(self,fname):

        self.crystalFile = fname
        self.particleCenters = []
        self.neighbs = {}


        # read in particle centers from csv
        with open(fname) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            first = True
            for row in reader:
                # first line of csv gives xmin, ymin, window dimensions, and bead radius
                if first:
                    self.beadRad = float(row[4])
                    self.xmin = float(row[0])
                    self.ymin = float(row[1])
                    self.windowSize = [float(row[2]),float(row[3])]
                    first = False
                    continue
                # after the first line, it's particle centers all the way down
                p = ( (float(row[0])-self.xmin) , 
                      (float(row[1])-self.ymin) )
                self.particleCenters.append(p)

        # read in neighbors
        i = self.crystalFile.find('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        neighbFile = "crysNeighbs/"+nameRoot+"_neighbs.csv"
        with open(neighbFile) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            for row in reader:
                # first element gives the particles whomst neighbors we will see
                partID = int(row[0])-1
                p = self.particleCenters[partID]
                
                # subsequent elements are that particle's neighbors
                myNeighbs = []
                for i in range(1,len(row)):
                    myNeighbs.append(self.particleCenters[int(row[i])-1])

                self.neighbs[p] = myNeighbs


    def polygonAt(self,p):
        # pick an N st the sidelength of an N-gon is ~= the length of a pixel
        N = int(np.ceil(2*np.pi*self.beadRad))
        delta = 2*np.pi/N
        return [ ( p[0]+self.beadRad*np.cos(n*delta), p[1]+self.beadRad*np.sin(n*delta) ) for n in range(N) ]
        
        
    def exclVolPolygonAt(self,p):
        #N = int(np.ceil(2*np.pi*self.beadRad))
        N = 64
        delta = 2*np.pi/N
        return [np.array([p[0]+2*self.beadRad*np.cos(n*delta) for n in range(N+1)]), np.array([p[1]+2*self.beadRad*np.sin(n*delta) for n in range(N+1)])]

    
    def showGrid(self):
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.scatter(*zip(*self.particleCenters),marker='.')
        plt.xlim(0,self.windowSize[0])
        plt.ylim(0,self.windowSize[1])
        #plt.xticks(np.arange(0, windowSize[0], step=1))
        #plt.yticks(np.arange(0, windowSize[1], step=1))
        #plt.grid(b=True,which='both',axis='both')
        plt.show()


    def showNeighbors(self,p):
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        for (x,y) in self.neighbs[p]:
            xvals = [x,p[0]]
            yvals = [y,p[1]]
            plt.plot(xvals,yvals,color='black')
        plt.scatter(*zip(*self.particleCenters),marker='.')
        plt.xlim(0,self.windowSize[0])
        plt.ylim(0,self.windowSize[1])
        plt.show()



    def freeSpace(self,particleID):
        # TODO return 0 for offgrid particles?

        center = self.particleCenters[particleID]
        nns = self.neighbs[center]

        crossingPts = [] #(x,y)
        crossingPairs = [] #(circle1,circle2)

        plt.plot(center[0],center[1],'*r')
        for nn in nns:
            plt.plot(nn[0],nn[1],'*r')

        # go over all pairs of circles
        for i in range(len(nns)):
            (x1,y1) = self.exclVolPolygonAt(nns[i])
            plt.plot(x1,y1)
            for j in range(i+1,len(nns)):
                (x2,y2) = self.exclVolPolygonAt(nns[j])
                
                myCrossingPts = myGeo.circIntersections(nns[i][0], nns[i][1], 2*self.beadRad, \
                                                        nns[j][0], nns[j][1], 2*self.beadRad)
                
                # no crossing points? skip
                if myCrossingPts == None:
                    continue
                
                # keep only the CLOSEST crossing point
                cross1 = (myCrossingPts[0],myCrossingPts[1])
                cross2 = (myCrossingPts[2],myCrossingPts[3])
                if dist(center,cross1) <= dist(center,cross2):
                    closestCrossingPt = cross1
                else:
                    closestCrossingPt = cross2
                
                # keep only crossing points that aren't contained in some other circle
                keepMe = True
                for k in range(len(nns)):
                    if k == i or k == j:
                        continue
                    if dist(nns[k],closestCrossingPt) < 2*self.beadRad:
                        keepMe = False
                        break

                if keepMe:
                    crossingPts.append(closestCrossingPt)
                    crossingPairs.append( (i,j) )
                    plt.plot(closestCrossingPt[0],closestCrossingPt[1],'*k')
        
        # sortKey is a function that returns a point's angle and distance relative to center
        # this will allow us to sort the crossing points in ccw order, which is a necessary pre-req for polyArea
        sortKey = myGeo.make_clockwiseangle_and_distance(center)
        myArea = myGeo.polyArea(sorted(crossingPts,key=sortKey))
        # go through each circle and cut out the appropriate segment
        # this is slightly inefficient but i think its ok
        for i in range(len(nns)):
            myPts = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPts.append(crossingPts[j])
            
            # i expect to always find 0 or 2 crossing points
            if len(myPts) == 0:
                continue
            if len(myPts) != 2:
                print("anna, a big problem has happened, please come fix it")
            
            vec1 = (myPts[0][0]-nns[i][0],myPts[0][1]-nns[i][1])
            vec2 = (myPts[1][0]-nns[i][0],myPts[1][1]-nns[i][1])
            dot = vec1[0]*vec2[0] + vec1[1]*vec2[1]
            cosine = dot/(4*self.beadRad*self.beadRad)
            theta = np.arccos(cosine)
            # segment area = (1/2) * (theta-sin(theta)) * R^2
            segArea = 0.5 * (theta-np.sin(theta)) * 4*self.beadRad*self.beadRad
            myArea = myArea - segArea

        plt.show()
        return myArea/(np.pi*self.beadRad*self.beadRad)


    def entropy(self,sbeadFile=None):
        tic = time.time()

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if sbeadFile == None:
            sbeadFile = nameRoot+'_Sbead.csv'


        S = 0 # total S
        numParts = 0 # number of particles counted
        buffer = self.beadRad*2

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            for i in range(len(self.particleCenters)):
                # don't count particles that are too close to the edge
                p = self.particleCenters[i]
                if p[0] < buffer or p[0] >= self.windowSize[0]-buffer or \
                p[1] < buffer or p[1] >= self.windowSize[1]-buffer:
                    continue
                
                numParts += 1

                [pID, freeArea] = self.freeSpace(i)
                Si = np.log(freeArea)
                S += Si
                writer.writerow([pID,Si])


        toc = time.time()
        return [S,numParts,toc-tic]

    def entropyParallel(self,numProc,sbeadFile=None):
        tic = time.time()

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
            
        if sbeadFile == None:
            sbeadFile = nameRoot+'_Sbead.csv'


        S = 0 # total S

        # pick out the particles to include in entropy calculation
        particlesInGrid = []
        buffer = self.beadRad*2
        for i in range(len(self.particleCenters)):
            # don't count particles that are too close to the edge
            p = self.particleCenters[i]
            if not(p[0] < buffer or p[0] >= self.windowSize[0]-buffer or \
            p[1] < buffer or p[1] >= self.windowSize[1]-buffer):
                particlesInGrid.append(i)
        numParts = len(particlesInGrid)

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            # get all snowflakes in parallel
            pool = mp.Pool(numProc)  
            pool_results = [pool.apply_async(self.freeSpace,args=[i]) for i in particlesInGrid]

            pool.close()
            pool.join()

            for r in pool_results:
                [pID,freeArea] = r.get()
                Si = np.log(freeArea)
                S += Si
                writer.writerow([pID,Si])

        toc = time.time()

        return [S,numParts,toc-tic]

   
    def volumeFraction(self):
        # TODO this would need to change for a hexagonal grid
        return len(self.occupiedPx) / self.windowSize[0] / self.windowSize[1]

   

def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

