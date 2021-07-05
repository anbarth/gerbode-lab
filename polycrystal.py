import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import multiprocessing as mp
import myGeo
import importlib
import os
from matplotlib import path
importlib.reload(myGeo)
import random

class Polycrystal:
    # class fields:
    #   xmin: minimum x value
    #   ymin: minimum y value TODO explain these better lol
    #   windowSize: a list that goes [width,height]. note that 0 is included and width and height are excluded.
    #   beadRad: bead radius in px
    #   particleCenters: 

    # constructor takes the input csv file
    # windowOverride=False (default) means count up entropy for all particles in polygon window (csv line 2)
    # windowOverride=True (useful for fake grainsplitting events) means countParticle is given line-by-line in the csv in the 3rd column
    def __init__(self,fname,neighbFile=None,windowOverride=False,radius=0):

        self.crystalFile = fname
        self.particleCenters = []
        self.neighbs = {}
        self.windowVertices = []
        # left, right, bot, top
        self.windowDims = [10000,-10000,10000,-1000] # dummy values to start
        # countParticle is a list of booleans, ith value says if particle i is in the window
        self.countParticle = []

        # read in particle centers from csv
        with open(fname) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            lineNum = 1
            for row in reader:
                # first line of csv gives bead radius
                if lineNum == 1:
                    self.beadRad = float(row[0])
                    if radius != 0:
                        self.beadRad = radius
                    lineNum += 1
                    continue
                # second line gives polygon window
                if lineNum == 2:                  
                    if len(row) % 2 != 0:
                        print("warning: row 2 should have an even number of entries, it doesn't")
                    for i in range(int(len(row)/2)):
                        myX = float(row[2*i])
                        myY = float(row[2*i+1])
                        self.windowVertices.append((myX,myY))

                        # adjust window bounds:
                        # left
                        if myX < self.windowDims[0]:
                            self.windowDims[0] = myX
                        # right
                        if myX > self.windowDims[1]:
                            self.windowDims[1] = myX
                        # bot
                        if myY < self.windowDims[2]:
                            self.windowDims[2] = myY
                        # top
                        if myY > self.windowDims[3]:
                            self.windowDims[3] = myY
                    lineNum += 1
                    continue
                        

                # after the first line, it's particle centers all the way down
                # col 0 gives particle ID, but that's just so that the csv is human-readable.
                # it's redundant information bc the particles need to be in particle ID order anyway

                # col 1 and 2 are (x,y) position
                p = ( (float(row[1])) , 
                      (float(row[2])) )
                self.particleCenters.append(p)
                
                # col 3 optionally says if this particle should be counted when calculating S
                # if this option is turned off, i'll just use all particles in the polygon window
                if windowOverride == True:
                    self.countParticle.append(bool(int(row[3])))

        # decide which beads are in the window and outside the window
        if windowOverride == False:
            windowPath = path.Path(self.windowVertices)
            self.countParticle = windowPath.contains_points(self.particleCenters)

        # if windowOverride is a string rather than a bool, interpret it as a file that specifies who's in and who's out
        if isinstance(windowOverride,str):
            with open(windowOverride) as windowFile:
                wreader = csv.reader(windowFile, delimiter=',')
                for wrow in wreader:
                    self.countParticle.append(bool(int(wrow[1])))


        # read in neighbors
        if neighbFile == None:
            i = self.crystalFile.find('/') # chop off any part of the name before a slash
            nameRoot = self.crystalFile[i+1:-4]
            neighbFile = "crysNeighbs/"+nameRoot+"_neighbs.csv"
        else:
            neighbFile = "crysNeighbs/"+neighbFile

        with open(neighbFile) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            for row in reader:
   
                # first element gives the particles whomst neighbors we will see
                partID = int(row[0])
                
                # subsequent elements are that particle's neighbors
                myNeighbs = []
                for i in range(1,len(row)):
                    myNeighbs.append(int(row[i]))

                # sortKey is a function that returns a point's angle and distance relative to p
                # this will allow us to sort the neighbors in ccw order, which is convenient in freeSpace
                sortKey = self.make_clockwiseangle_and_distance_byID(partID)
                self.neighbs[partID] = sorted(myNeighbs,key=sortKey)
                
    
    def foo(self):
        for i in range(len(self.particleCenters)):
            for j in range(i+1,len(self.particleCenters)):
                myR = dist(self.particleCenters[i],self.particleCenters[j])/2
                if myR <= 11:
                    print(i+1,j+1,myR)


    def show(self):
        fig, ax = plt.subplots()
        ax.set_aspect(1)

        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]
            circ = plt.Circle(p, self.beadRad, facecolor='gray',edgecolor=None, alpha=0.3)
            ax.add_artist(circ)
    
            if self.countParticle[i]:
                circ = plt.Circle(p, self.beadRad/5, facecolor='k', edgecolor=None)
                ax.add_artist(circ)
            
            #plt.text(p[0]+self.beadRad/15,p[1]+self.beadRad/15,str(i+1))

        plt.xlim(self.windowDims[0]-2*self.beadRad,self.windowDims[1]+2*self.beadRad)
        plt.ylim(self.windowDims[2]-2*self.beadRad,self.windowDims[3]+2*self.beadRad)

        plt.show()
        #strPos = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        #nameRoot = self.crystalFile[strPos+1:-4]
        #fig.savefig(nameRoot+"_img.png",dpi=900)
        plt.close(fig)
        


    def showNeighbors(self,pID):

        fig, ax = plt.subplots()
        ax.set_aspect(1)
        
        for i in range(len(self.particleCenters)):
            q = self.particleCenters[i]
            circ = plt.Circle(q, self.beadRad, facecolor='gray',edgecolor=None, alpha=0.3)
            ax.add_artist(circ)
    
            if self.countParticle[i]:
                circ = plt.Circle(q, self.beadRad/20, facecolor='k', edgecolor=None)
                ax.add_artist(circ)

        # keep track of neighbors' positions just so i know how big to make the window
        x_list = []
        y_list = []    
        for nnID in self.neighbs[pID]:
            (x,y) = self.particleCenters[nnID-1]
            x_list.append(x)
            y_list.append(y)
            circ = plt.Circle((x,y), 2*self.beadRad, facecolor=(1, 0, 0, 0.25),edgecolor='r')
            ax.add_artist(circ)
            plt.text(x+self.beadRad/15,y+self.beadRad/15,str(nnID))

        myLeft = min(x_list)
        myRight = max(x_list)
        myBot = min(y_list)
        myTop = max(y_list)

        plt.xlim(myLeft-self.beadRad,myRight+self.beadRad)
        plt.ylim(myBot-self.beadRad,myTop+self.beadRad)
        plt.show()
        plt.close(fig)

    
    # returns the area and shape of the particle specified by particleID
    def freeSpace(self,particleID,imgMaking=False):
        thisWentSmoothly = True

        # return 0 for particles that shouldn't be counted
        if not self.countParticle[particleID-1]:
            return 0

        center = self.particleCenters[particleID-1]
        nnIDs = self.neighbs[particleID]
        nns = [self.particleCenters[nnID-1] for nnID in nnIDs]

        # (x,y) for crossing points between neighbors' excl area circles
        crossingPts = [] 
        # (circle1,circle2) specifying, for each entry in crossingPts, which neighbors' circles are crossing
        crossingPairs = [] 

        # go over all pairs of circles
        for i in range(len(nns)):
            for j in range(i+1,len(nns)): 
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
        
        # check for and correct nearest neighbor network issues
        for i in range(len(nns)):
            # find the crossing points on this circle
            myPtInds = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPtInds.append(j)

            # there's an issue if there's >2 crossing points
            if len(myPtInds) <= 2:
                continue
            
            if thisWentSmoothly:
                print(particleID,"has extraneous neighbors")
                thisWentSmoothly = False
            #print("============")
            #print(nnIDs[i])
            # find the two crossing points closest to the particle's center
            crossingPtInd1 = -1
            crossingPtInd2 = -1
            crossingPtDist1 = np.inf
            crossingPtDist2 = np.inf
            for j in myPtInds:
                myDist = dist(crossingPts[j],center)
                if myDist < crossingPtDist1:
                    crossingPtInd2 = crossingPtInd1
                    crossingPtDist2 = crossingPtDist1
                    crossingPtInd1 = j
                    crossingPtDist1 = myDist
                elif myDist < crossingPtDist2:
                    crossingPtInd2 = j
                    crossingPtDist2 = myDist

            # now, get rid of all the crossing points except the closest two
            myPtInds = sorted(myPtInds,reverse=True) # modifying a list, so work back to front
            for j in myPtInds:
                if j != crossingPtInd1 and j != crossingPtInd2:
                    #print(crossingPts[j])
                    #print(nnIDs[crossingPairs[j][0]],nnIDs[crossingPairs[j][1]])
                    crossingPts = np.delete(crossingPts,j,0)
                    crossingPairs = np.delete(crossingPairs,j,0)

        # check for stragglers
        #print(crossingPairs)
        #print(crossingPts)
        for i in range(len(nns)):
            # find the crossing points on this circle
            myPtInds = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPtInds.append(j)

            if len(myPtInds) == 1:
                j = myPtInds[0]
                #print("straggler:")
                #print(crossingPts[j])
                #print(nnIDs[crossingPairs[j][0]],nnIDs[crossingPairs[j][1]])
                crossingPts = np.delete(crossingPts,j,0)
                crossingPairs = np.delete(crossingPairs,j,0)

        if len(crossingPts) == 0:
            print("free area appears to be 0 or negative for "+str(particleID)+", go check it out")
            return [particleID,0,[center[0]],[center[1]]]

        # find the area of the polygon bounded by the crossing points
        # first, find a point (x_inside, y_inside) that is inside the polygon
        x_inside = np.average([crossPt[0] for crossPt in crossingPts])
        y_inside = np.average([crossPt[1] for crossPt in crossingPts])
        # now define sortKey, a function that returns a point's angle and distance relative to (x_inside,y_inside)
        # this allows us to sort the crossing points in ccw order, which is a necessary pre-req for polyArea
        sortKey = myGeo.make_clockwiseangle_and_distance((x_inside,y_inside))
        myArea = myGeo.polyArea(sorted(crossingPts,key=sortKey))

        # find the segments of excluded area that protrude into the polygon
        # these lists will store the (x,y) points that define the free space's perimeter
        freeSpaceCurveX = []
        freeSpaceCurveY = []

        # go through each excluded area circle and cut out the appropriate segment
        for i in range(len(nns)):

            # find the crossing points on this circle
            myPts = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPts.append(crossingPts[j])
            
            # i expect to always find 0 or 2 crossing points
            if len(myPts) == 0: 
                continue
            if len(myPts) != 2 or len(myPts) == 1:
                print(particleID,"somehow still has an issue, come fix it")
                thisWentSmoothly = False

            # vectors pointing from this neighbor to each crossing pt
            vec1 = (myPts[0][0]-nns[i][0],myPts[0][1]-nns[i][1])
            vec2 = (myPts[1][0]-nns[i][0],myPts[1][1]-nns[i][1])

            # angle of each vector, relative to the positive x axis
            theta1 = np.arctan2(vec1[1],vec1[0])
            theta2 = np.arctan2(vec2[1],vec2[0])

            # arctan2's range is [-pi,pi], meaning there's an awkward jump at pi
            # adjust theta1, theta2 so that they do not straddle the jump
            # i jump through all these hoops so that later on, you can plot a 
            # pretty, continuous range of thetas between theta1 and theta2
            thetann = np.arctan2(nns[i][1]-center[1],nns[i][0]-center[0])
            if theta1 >= thetann and theta1 <= np.pi:
                theta1 = theta1 - 2*np.pi
            if theta2 >= thetann and theta2 <= np.pi:
                theta2 = theta2 - 2*np.pi
            if theta1 >= -np.pi and theta1 < -2*np.pi+thetann:
                theta1 = theta1 + 2*np.pi
            if theta2 >= -np.pi and theta2 < -2*np.pi+thetann:
                theta2 = theta2 + 2*np.pi

            # positive angle between vec1 and vec2
            theta = abs(theta1-theta2)
            # segment area = (1/2) * (theta-sin(theta)) * R^2
            segArea = 0.5 * (theta-np.sin(theta)) * 4*self.beadRad*self.beadRad
            myArea = myArea - segArea

            # add segment points to freeSpaceCurve
            thetaMin = min(theta1,theta2)
            thetaMax = max(theta1,theta2)
            numSteps = 50 # idk just some number i like

            freeSpaceCurveX.extend([nns[i][0]+2*self.beadRad*np.cos(t) for t in np.linspace(thetaMin,thetaMax,numSteps)])
            freeSpaceCurveY.extend([nns[i][1]+2*self.beadRad*np.sin(t) for t in np.linspace(thetaMin,thetaMax,numSteps)])

        freeArea = myArea/(np.pi*self.beadRad*self.beadRad)
        if freeArea > 0.5:
            print(particleID,"has rather large free area, worth checking out")
            thisWentSmoothly = False

        if not(thisWentSmoothly) and imgMaking:
            fig, ax = plt.subplots()
            ax.set_aspect(1)

            circ = plt.Circle(center, self.beadRad/60, color='black')
            ax.add_artist(circ)
            plt.text(center[0]+self.beadRad/15,center[1]+self.beadRad/15,str(particleID))

            for i in range(len(nns)):
                circ = plt.Circle(nns[i], self.beadRad/60, color='black')
                ax.add_artist(circ)
                plt.text(nns[i][0]+self.beadRad/15,nns[i][1]+self.beadRad/15,str(nnIDs[i]))
                circ = plt.Circle(nns[i], 2*self.beadRad, color='gray', alpha=0.3)
                ax.add_artist(circ)

            plt.plot(freeSpaceCurveX,freeSpaceCurveY)
            plt.xlim(min(freeSpaceCurveX)-self.beadRad/3,max(freeSpaceCurveX)+self.beadRad/3)
            plt.ylim(min(freeSpaceCurveY)-self.beadRad/3,max(freeSpaceCurveY)+self.beadRad/3)

            strPos = self.crystalFile.find('/') # chop off any part of the name before a slash
            nameRoot = self.crystalFile[strPos+1:-4]
            myFileName = "badFreeSpace/"+nameRoot+"_freespace_"+str(particleID)+".png"
            try:
                fig.savefig(myFileName, dpi=900)
            except FileNotFoundError:
                strPos = nameRoot.rfind('/')
                os.makedirs("badFreeSpace/"+nameRoot[0:strPos])
                fig.savefig(myFileName, dpi=900)
            plt.close(fig)

        return [particleID,myArea/(np.pi*self.beadRad*self.beadRad),freeSpaceCurveX,freeSpaceCurveY]


    def drawFreeSpace(self,particleID,show=True,save=False):
        p = self.particleCenters[particleID-1]
        (pID,freeArea,freeSpaceCurveX,freeSpaceCurveY) = self.freeSpace(particleID)

        fig, ax = plt.subplots()
        ax.set_aspect(1)

        circ = plt.Circle(p, self.beadRad/60, color='black')
        ax.add_artist(circ)

        nnIDs = self.neighbs[particleID]
        nns = [self.particleCenters[nnID-1] for nnID in nnIDs]
        for i in range(len(nns)):
            circ = plt.Circle(nns[i], self.beadRad/60, color='black')
            ax.add_artist(circ)
            plt.text(nns[i][0]+self.beadRad/15,nns[i][1]+self.beadRad/15,str(nnIDs[i]))
            circ = plt.Circle(nns[i], 2*self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.plot(freeSpaceCurveX,freeSpaceCurveY)
        plt.xlim(min(freeSpaceCurveX)-self.beadRad/3,max(freeSpaceCurveX)+self.beadRad/3)
        plt.ylim(min(freeSpaceCurveY)-self.beadRad/3,max(freeSpaceCurveY)+self.beadRad/3)
        plt.show()
        plt.close(fig)
        return freeArea

    
    def entropy(self,sbeadFile=None,imgFile=None):
        tic = time.time()
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        plt.xlim(self.windowDims[0]-2*self.beadRad,self.windowDims[1]+2*self.beadRad)
        plt.ylim(self.windowDims[2]-2*self.beadRad,self.windowDims[3]+2*self.beadRad)

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if sbeadFile == None:
            sbeadFile = 'freeSpaceCSVs/'+nameRoot+'_freeSpaces.csv'
        if imgFile == None:
            imgFile = 'snowflakes/'+nameRoot+'_snowflakes.png'
        cmap = cm.get_cmap('viridis')

        S = 0 # total S
        numParts = 0 # number of particles counted
        buffer = self.beadRad*2

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            for i in range(len(self.particleCenters)):
                p = self.particleCenters[i]

                # draw
                circ = plt.Circle(p, self.beadRad, facecolor='#b8b8b8',edgecolor='black',linewidth=0,alpha=1,zorder=0)
                ax.add_artist(circ)

                # don't count particles that are outside the window
                if not self.countParticle[i]:
                    continue
                numParts += 1

                # draw a black dot to indicate this particle is included
                circ = plt.Circle(p, self.beadRad/30, facecolor='k', edgecolor=None,zorder=10)
                ax.add_artist(circ) 

                # find the free area -- note that particleID = i+1 !!
                (pID,freeArea,freeSpaceCurveX,freeSpaceCurveY) = self.freeSpace(i+1,imgMaking=True)

                if freeArea <= 0:
                    Si = 0
                else:
                    Si = np.log(freeArea)

                S += Si

                writer.writerow([pID,freeArea]) # record in Sbead file

                # draw free area
                # cmap takes a number in [0,1) to a color
                rgb = cmap((freeArea)*25)[0:3]
                #ax.fill(freeSpaceCurveX,freeSpaceCurveY, facecolor=rgb,edgecolor='black',lw=0.15)

                # to plot free space at 3x its size, use this code instead
                cen = centroid(freeSpaceCurveX,freeSpaceCurveY)
                biggerCurveX = 2*(np.array(freeSpaceCurveX)-cen[0])+cen[0]
                biggerCurveY = 2*(np.array(freeSpaceCurveY)-cen[1])+cen[1]
                ax.fill(biggerCurveX,biggerCurveY, facecolor=rgb,edgecolor='black',lw=0.15,zorder=5)
        

        fig.savefig(imgFile, dpi=900)
        plt.close(fig)
        toc = time.time()
        return [S,numParts,toc-tic]

    def entropyMC(self,numTrials):
        successfulTrials = 0
        for i in range(numTrials):
            successfulTrials += self.MCtrial()
        return successfulTrials/numTrials

    def MCtrial(self):
        trialCenters = []
        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]

            # x0 and y0 are in [-2beadRad,+2beadRad)
            x0 = random.random()*4*self.beadRad - 2*self.beadRad
            y0 = random.random()*4*self.beadRad - 2*self.beadRad
            
            trialCenters.append( (p[0]+x0, p[1]+y0) )
        
        for pID in range(1,len(trialCenters)+1):
            myNeighbs = self.neighbs[pID]
            for nnID in myNeighbs:
                if dist(trialCenters[pID-1],trialCenters[nnID-1]) < self.beadRad*2:
                    return 0
        
        return 1


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
            p = self.particleCenters[i]

            # don't count particles that are outside the window
            if self.countParticle[i]:
                particlesInGrid.append(i+1)

        numParts = len(particlesInGrid)

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            # get all snowflakes in parallel
            pool = mp.Pool(numProc)  
            pool_results = [pool.apply_async(self.freeSpace,args=[pID]) for pID in particlesInGrid]

            pool.close()
            pool.join()

            for r in pool_results:
                [pID,freeArea] = r.get()
                Si = np.log(freeArea)
                S += Si
                writer.writerow([pID,Si])

        toc = time.time()

        return [S,numParts,toc-tic]

    def areaFraction(self):  
        # area occupied = N pi r^2
        # TODO this isn't accounting for particles that fall partially outside the window
        windowPath = path.Path(self.windowVertices)
        inWindow = windowPath.contains_points(self.particleCenters)
        numParts = sum(inWindow)
        areaOccupied = numParts*np.pi*self.beadRad*self.beadRad
        totalArea = myGeo.polyArea(self.windowVertices)
        return areaOccupied/totalArea
    
    # im so sorry this god-awful function exists twice in my code base
    # it's once here and there's also an extremely similar function in myGeo
    # the myGeo one is older and doesn't really need to exist if this one's here
    # i would rather this ugly thing be in myGeo, where i don't have to look at it
    # but i can't figure out how to sort neighbors in ccw without this here
    def make_clockwiseangle_and_distance_byID(self,originID):
        origin = self.particleCenters[originID-1]
        def clockwiseangle_and_distance_byID(pID):
            refvec = [0,1]
            # Vector between point and the origin: v = p - o
            point = self.particleCenters[pID-1]
            vector = [point[0]-origin[0], point[1]-origin[1]]
            # Length of vector: ||v||
            lenvector = np.hypot(vector[0], vector[1])
            # If length is zero there is no angle
            if lenvector == 0:
                return (-np.pi, 0)
            # Normalize vector: v/||v||
            normalized = [vector[0]/lenvector, vector[1]/lenvector]
            dotprod  = normalized[0]*refvec[0] + normalized[1]*refvec[1]     # x1*x2 + y1*y2
            diffprod = refvec[1]*normalized[0] - refvec[0]*normalized[1]     # x1*y2 - y1*x2
            angle = np.arctan2(diffprod, dotprod)
            # Negative angles represent counter-clockwise angles so we need to subtract them 
            # from 2*pi (360 degrees)
            if angle < 0:
                return 2*np.pi+angle, lenvector
            # I return first the angle because that's the primary sorting criterium
            # but if two vectors have the same angle then the shorter distance should come first.
            return angle, lenvector
        return clockwiseangle_and_distance_byID

def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

def centroid(x_coords,y_coords):
    _len = len(x_coords)
    centroid_x = sum(x_coords)/_len
    centroid_y = sum(y_coords)/_len
    return [centroid_x, centroid_y]