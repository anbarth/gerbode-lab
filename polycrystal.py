import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
import multiprocessing as mp
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
                if windowOverride:
                    self.countParticle.append(bool(int(row[3])))

        # decide which beads are in the window and outside the window
        if not windowOverride:
            windowPath = path.Path(self.windowVertices)
            self.countParticle = windowPath.contains_points(self.particleCenters)

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
                # TODO its super confusing that matlab IDs are pID+1 :(
                # first element gives the particles whomst neighbors we will see
                partID = int(row[0])
                p = self.particleCenters[partID-1]
                
                # subsequent elements are that particle's neighbors
                myNeighbs = []
                for i in range(1,len(row)):
                    myNeighbs.append(self.particleCenters[int(row[i])-1])

                # sortKey is a function that returns a point's angle and distance relative to p
                # this will allow us to sort the neighbors in ccw order, which is convenient in freeSpace
                sortKey = myGeo.make_clockwiseangle_and_distance(p)
                self.neighbs[p] = sorted(myNeighbs,key=sortKey)
    
    def show(self):
        fig, ax = plt.subplots()
        ax.set_aspect(1)

        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]
            circ = plt.Circle(p, self.beadRad, facecolor='gray',edgecolor=None, alpha=0.3)
            ax.add_artist(circ)
    
            if self.countParticle[i]:
                circ = plt.Circle(p, self.beadRad/10, facecolor='k', edgecolor=None)
                ax.add_artist(circ)
                #plt.text(p[0]+self.beadRad/15,p[1]+self.beadRad/15,str(i+1))

        plt.xlim(self.windowDims[0],self.windowDims[1])
        plt.ylim(self.windowDims[2],self.windowDims[3])


        plt.show()


    def showNeighbors(self,pID):
        p = self.particleCenters[pID-1]
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        
        for i in range(len(self.particleCenters)):
            q = self.particleCenters[i]
            circ = plt.Circle(q, self.beadRad, facecolor='gray',edgecolor=None, alpha=0.3)
            ax.add_artist(circ)
    
            if self.countParticle[i]:
                circ = plt.Circle(q, self.beadRad/20, facecolor='k', edgecolor=None)
                ax.add_artist(circ)
            
            plt.text(q[0]+self.beadRad/15,q[1]+self.beadRad/15,str(i+1))

        for (x,y) in self.neighbs[p]:
            circ = plt.Circle((x,y), 2*self.beadRad, facecolor=(0, 0, 0, 0),edgecolor='r')
            ax.add_artist(circ)

        myLeft = min([x for (x,y) in self.neighbs[p]])
        myRight = max([x for (x,y) in self.neighbs[p]])
        myBot = min([y for (x,y) in self.neighbs[p]])
        myTop = max([y for (x,y) in self.neighbs[p]])

        plt.xlim(myLeft-self.beadRad,myRight+self.beadRad)
        plt.ylim(myBot-self.beadRad,myTop+self.beadRad)
        plt.show()

    
    # returns the area and shape of the particle specified by particleID
    def freeSpace(self,particleID):
        # return 0 for particles that shouldn't be counted
        if not self.countParticle[particleID-1]:
            return 0

        center = self.particleCenters[particleID-1]
        nns = self.neighbs[center]

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
        
        # find the area of the polygon bounded by the crossing points
        # sortKey is a function that returns a point's angle and distance relative to center
        # this allows us to sort the crossing points in ccw order, which is a necessary pre-req for polyArea
        sortKey = myGeo.make_clockwiseangle_and_distance(center)
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
            if len(myPts) != 2:
                print(particleID,"has extraneous neighbors, come fix it")
            
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
            numSteps = int( (thetaMax-thetaMin)/(np.pi/180) ) # number of steps st every step covers 1 degree

            freeSpaceCurveX.extend([nns[i][0]+2*self.beadRad*np.cos(t) for t in np.linspace(thetaMin,thetaMax,numSteps)])
            freeSpaceCurveY.extend([nns[i][1]+2*self.beadRad*np.sin(t) for t in np.linspace(thetaMin,thetaMax,numSteps)])

        freeArea = myArea/(np.pi*self.beadRad*self.beadRad)
        if freeArea > 0.5:
            print(particleID,"has rather large free area, worth checking out")

        return [particleID,myArea/(np.pi*self.beadRad*self.beadRad),freeSpaceCurveX,freeSpaceCurveY]


    def drawFreeSpace(self,particleID):
        p = self.particleCenters[particleID-1]
        (pID,freeArea,freeSpaceCurveX,freeSpaceCurveY) = self.freeSpace(particleID)

        fig, ax = plt.subplots()
        ax.set_aspect(1)

        plt.scatter(p[0],p[1])
        #circ = plt.Circle(p, self.beadRad, color='gray', alpha=0.3)
        #ax.add_artist(circ)

        for i in range(len(self.particleCenters)):
            q = self.particleCenters[i]
            plt.text(q[0]+self.beadRad/15,q[1]+self.beadRad/15,str(i+1))

        for nn in self.neighbs[p]:
            plt.scatter(nn[0],nn[1])
            circ = plt.Circle(nn, 2*self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.plot(freeSpaceCurveX,freeSpaceCurveY)
        plt.show()
        return freeArea


    def entropy(self,sbeadFile=None,imgFile=None):
        tic = time.time()
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        #plt.xlim(140,460)
        #plt.ylim(140,460)

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if sbeadFile == None:
            sbeadFile = nameRoot+'_freeSpaces.csv'
        if imgFile == None:
            imgFile = nameRoot+'_snowflakes.png'
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
                circ = plt.Circle(p, self.beadRad/15, facecolor='k', edgecolor=None,zorder=10)
                ax.add_artist(circ) 

                # find the free area -- note that particleID = i+1 !!
                (pID,freeArea,freeSpaceCurveX,freeSpaceCurveY) = self.freeSpace(i+1)

                if freeArea == 0:
                    Si = 0
                    print("free area appears to be 0 for "+str(i+1)+", go check it out")
                else:
                    Si = np.log(freeArea)

                S += Si

                writer.writerow([pID,freeArea]) # record in Sbead file

                # draw free area
                # pick a color (freeArea = 0 --> blue; freeArea >= pi R^2/6 --> yellow)
                # cmap takes a number 0 to 1
                rgb = cmap((freeArea-0.035)*6.9)[0:3]
                #print(freeArea)
                #cen = centroid(freeSpaceCurveX,freeSpaceCurveY)
                #biggerCurveX = 3*(np.array(freeSpaceCurveX)-cen[0])+cen[0]
                #biggerCurveY = 3*(np.array(freeSpaceCurveY)-cen[1])+cen[1]
                ax.fill(freeSpaceCurveX,freeSpaceCurveY, facecolor=rgb,edgecolor='black',lw=0.15)
                #ax.fill(biggerCurveX,biggerCurveY, facecolor=rgb,edgecolor='black',lw=0.5,zorder=5)
        

        fig.savefig(imgFile, dpi=900)
        plt.close(fig)
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
    

def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

def centroid(x_coords,y_coords):
    _len = len(x_coords)
    centroid_x = sum(x_coords)/_len
    centroid_y = sum(y_coords)/_len
    return [centroid_x, centroid_y]