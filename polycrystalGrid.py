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

####   polycrystalgrid is the polycrystal class's dark past...
####   it's like polycrystal, except all positions exist on a grid.
####   and so then you calculate free area by literally counting up pixels
####   like some kind of caveman.
####   the comments are incomplete and i'm not gonna fix them because you should just 
####   use polycrystal anyway. this class isnt useful for anything
####   and i am including it more for the memories than anything else. 
####   anna barth 2021

class PolycrystalGrid:
    # class fields:
    #   xmin: minimum x value
    #   ymin: minimum y value TODO explain these better lol
    #   gridSize: a list that goes [width,height]. note that 0 is included and width and height are excluded.
    #   beadRad: bead radius in px
    #   particleCenters: 
    #   occupiedPx: a list of grid spots occupied by (x,y)s
    #   resolution: 
    #   psi6: does the file contain psi6 magnitude data?

    # constructor takes the input csv file
    def __init__(self,fname,rad=0,usePsi6=False,useNeighbs=False):

        self.crystalFile = fname
        self.particleCenters = []
        self.psi6dict = {}
        self.occupiedPx = []
        self.exclShape = []
        self.beadShape = []
        self.snowflakeShape = []
        self.usePsi6 = usePsi6
        self.psi6cutoff = 0.98
        self.useNeighbs = useNeighbs
        self.neighbs = {}


        # read in particle centers (and optionally, |psi6|) from csv
        with open(fname) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            first = True
            for row in reader:
                # TODO resolution stuff would probably need to change for a hexagonal grid
                # first line of csv gives xmin, ymin, grid dimensions, and bead radius
                if first:
                    resolution = 0
                    if rad == 0:
                        resolution = 1
                        self.beadRad = round( float(row[4])*resolution )
                    else:
                        resolution = rad/float(row[4])
                        self.beadRad = rad
                    self.xmin = int(row[0])
                    self.ymin = int(row[1])
                    self.gridSize = [round( int(row[2])*resolution ),round( int(row[3])*resolution )]
                    
                    first = False
                    continue
                # after the first line, it's particle centers all the way down
                p = (round( (float(row[0])-self.xmin)*resolution ), 
                     round( (float(row[1])-self.ymin)*resolution ))
                self.particleCenters.append(p)
                if self.usePsi6:
                    self.psi6dict[p] = float(row[-1])
                
        # particles and exclused volumes always have the same shapes
        # so set those shapes once at the beginning and recycle them
        self.setBeadShape() # do this one first, bc setExclShape relies on setBeadShape
        self.setExclShape() 

        # fill in all the occupied pixels
        for (x0,y0) in self.particleCenters:
            self.occupiedPx.extend(self.pxOccupiedByParticle((x0,y0)))
        self.occupiedPx = sorted(self.occupiedPx)

        if self.usePsi6:
            self.setSnowflakeShape()

        if self.useNeighbs:
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

    # sets the shape of a particle centered at (0,0), so you can translate it to any other position
    # ^ that's what pxOccupiedByParticle does
    # WARNING i am assuming your window is large enough to accomodate a particle at its center
    # if it's not, idk do something else
    def setBeadShape(self):
        # suppose we put a particle at the center of the grid, 
        # where it's (hopefully) far from any annoying edges that might obscure its beautiful round shape
        x0 = int(self.gridSize[0]/2)
        y0 = int(self.gridSize[1]/2)

        # find px included in this bead by starting at the center and spiraling out until you're outside the bead radius
        pxToCheck = [(x0,y0)]
        occupied = []
        for (x,y) in pxToCheck:
            # TODO this should be either inclusive or exclusive idk
            # been mostly using inclusive but i question myself :/
            if dist((x0,y0),(x,y)) <= self.beadRad:
                occupied.append((x-x0,y-y0))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)

        self.beadShape = occupied
        return occupied


    # sets the shape of the area excluded by a particle centered at (0,0), so you can translate it to any other position
    # ^ that's what pxExcludedByParticle does
    # WARNING i am assuming your window is large enough to accomodate a particle at its center
    # if it's not, idk do something else
    # also note that you need to run setBeadShape first!
    def setExclShape(self):
        # suppose we put a particle at the center of the grid, 
        # where it's (hopefully) far from any annoying edges that might obscure its beautiful round shape
        x0 = int(self.gridSize[0]/2)
        y0 = int(self.gridSize[1]/2)
        
        # get all the spots occupied by a bead at (x0,y0)
        # sort that shit bc we're gonna be searching it a LOT
        pxList = sorted(self.pxOccupiedByParticle((x0,y0)))

        # find px excluded by this bead by starting at the center and spiraling out
        pxToCheck = [(x0,y0)]
        excluded = []
        for (x,y) in pxToCheck:
            if self.conflict((x,y),pxList):
                excluded.append((x-x0,y-y0))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)

        self.exclShape = excluded
        return excluded

    # neatParticlePos is the position of a particle with perfect or near-perfect |psi6|=1
    # in the future i'd like to be able to just look up a particle like this automatically, ie look for max psi6
    # or the average |psi6| among particles below the cutoff
    def setSnowflakeShape(self):
        # when its done, try running testSimple and get entropy for annasCoolTest10.csv
        for i in range(len(self.particleCenters)):
            buffer = self.beadRad*2
            p = self.particleCenters[i]
            if p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
               p[1] < buffer or p[1] >= self.gridSize[1]-buffer:
                continue
            if self.psi6dict[p] >= self.psi6cutoff:
                [pID,freePx] = self.freeSpace(i)
                snowflake = [(x-p[0],y-p[1]) for (x,y) in freePx]
                self.snowflakeShape = snowflake
                return snowflake
        print("failed to find a special snowflake </3")
        return [] # failure :(


    def getNeighbors(self,p):
        # TODO this would need to change for a hexagonal grid
        (x,y) = p
        neighbors = []
        # remember: we include 0 and exclude gridSize
        if x > 0:
            neighbors.append((x-1,y))
        if y > 0:
            neighbors.append((x,y-1))
        if x < self.gridSize[0]-1:
            neighbors.append((x+1,y))
        if y < self.gridSize[1]-1:
            neighbors.append((x,y+1))
        return neighbors
    
    def showGrid(self):
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.scatter(*zip(*self.occupiedPx),marker='.')
        plt.xlim(0,self.gridSize[0])
        plt.ylim(0,self.gridSize[1])
        #plt.xticks(np.arange(0, gridSize[0], step=1))
        #plt.yticks(np.arange(0, gridSize[1], step=1))
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
        plt.scatter(*zip(*self.occupiedPx),marker='.')
        plt.xlim(0,self.gridSize[0])
        plt.ylim(0,self.gridSize[1])
        plt.show()

    def showGridNew(self):
        scale = 5
        # initialize a white grid
        imgArr = np.ones((self.gridSize[1]*scale,self.gridSize[0]*scale,3),dtype=np.uint8) * 255
        # set all the occupied px to blue
        for (x,y) in self.occupiedPx:
            for i in range(scale):
                for j in range(scale):
                    imgArr[y*scale+j,x*scale+i] = [24,96,148]
        # TODO you should write the img to a csv or smtg
        img = PIL.Image.fromarray(imgArr,'RGB')
        img.show()


    def pxExcludedByParticle(self,p):
        pxExcl = []
        for (x,y) in self.exclShape:
            xpos = x+p[0]
            ypos = y+p[1]
            if xpos >= 0 and xpos < self.gridSize[0] and ypos >= 0 and ypos < self.gridSize[1]:
                pxExcl.append((xpos,ypos))
        return pxExcl

    def pxOccupiedByParticle(self,p):
        pxOcc = []
        for (x,y) in self.beadShape:
            xpos = x+p[0]
            ypos = y+p[1]
            if xpos >= 0 and xpos < self.gridSize[0] and ypos >= 0 and ypos < self.gridSize[1]:
                pxOcc.append((xpos,ypos))
        return pxOcc


    # p=(x0,y0) is a test position. does it conflict with the positions in pxList?
    # TODO offgrid positions should never be available
    def conflict(self,p,pxList):
        # special case: if pxList empty, no conflict
        if len(pxList) == 0:
            return False

        # pretend to put a bead at p
        # TODO can probly come up w a better name lol
        dangerZone = self.pxOccupiedByParticle(p)

        # now see if any of the positions occupied by the bead are in pxList
        for (x,y) in dangerZone:
            index = bisect.bisect_left(pxList,(x,y))
            if index != len(pxList) and pxList[index] == (x,y): # fancy version of (x,y) in pxList
                return True
        return False
    
    def isAvailable(self,p):
        return not self.conflict(p,self.occupiedPx)

    # TODO offgrid positions should never be available
    # TODO it'd be more elegant to use self.conflict, but also maybe slower
    # same as isAvailable, but don't include the positions listed in ignoreMe
    def isAvailableIgnore(self,p,ignoreMe):
        dangerZone = self.pxOccupiedByParticle(p)

        for (x,y) in dangerZone:
            if (x,y) not in ignoreMe:
                index = bisect.bisect_left(self.occupiedPx,(x,y))
                if index != len(self.occupiedPx) and self.occupiedPx[index] == (x,y):
                    return False
        return True

    def freeSpace(self,particleID):
        particleCenter = self.particleCenters[particleID]
        particle = self.pxOccupiedByParticle(particleCenter) # particle includes all (x,y) to ignore
        freePx = []
        pxToCheck = [particleCenter]
        for (x,y) in pxToCheck:
            if self.isAvailableIgnore((x,y),particle):
                freePx.append((x,y))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)
        return [particleID,freePx]

    def freeSpaceGeo(self,particleID):
        # TODO return 0 for offgrid particles?

        center = self.particleCenters[particleID]
        print("center",center)
        nns = self.neighbs[center]

        crossingPts = [] #(x,y)
        crossingPairs = [] #(circle1,circle2)

        plt.plot(center[0],center[1],'*r')
        for nn in nns:
            plt.plot(nn[0],nn[1],'*r')

        # go over all pairs of circles
        for i in range(len(nns)):
            print("neighbor",i,nns[i])
            (x1,y1) = self.exclVolPolygonAt(nns[i])
            plt.plot(x1,y1)
            for j in range(i+1,len(nns)):
                #(x2,y2) = self.exclVolPolygonAt(nns[j])
                
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
        

        myArea = myGeo.polyArea(crossingPts)
        print("polygon area",myArea)
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
            print("segment for neighbor ",i,"area",segArea)
            myArea = myArea - segArea
        
        #print(myArea)
        plt.show()
        return [particleID, myArea/(np.pi*self.beadRad*self.beadRad)]


    def isAvailPoly(self,oldCenter,newCenter):
        # particle can always exist in the spot it originally is
        # we need to specify this bc the crystal might actually have some stupid little bits over overlap
        if oldCenter == newCenter:
            return True
        nns = self.neighbs[oldCenter]

        # loops over nearest neighbors
        for nnCenter in nns:
            d = dist(newCenter,nnCenter)
            if d < 2*self.beadRad:
                return False
        # made it through all the neighbors? you're good to go
        return True


    def freeSpacePoly(self,particleID):
        particleCenter = self.particleCenters[particleID]
        pxToCheck = [particleCenter]
        freePx = []
        for (x,y) in pxToCheck:
            if self.isAvailPoly(particleCenter,(x,y)):
                freePx.append((x,y))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)
        return [particleID,freePx]

    def freeAreaPoly(self,particleID):
        (pID,freePx) = self.freeSpacePoly(particleID)
        return len(freePx)/len(self.beadShape)

    def freeSpaceMC(self,particleID):
        particleCenter = self.particleCenters[particleID]

        # randomly guess points within 1 diameter = 2 beadRad of particleCenter
        maxTrials = 500
        numTrials = 0
        freePts = []
        while numTrials < maxTrials:
            # x0 and y0 are in [-2beadRad,+2beadRad)
            #x0 = random.random()*4*self.beadRad - 2*self.beadRad
            #y0 = random.random()*4*self.beadRad - 2*self.beadRad
            x0 = random.random()*2*self.beadRad - self.beadRad
            y0 = random.random()*2*self.beadRad - self.beadRad
            # reject points that are not in the circle (excluding boundary but idt it matters)
            if x0*x0 + y0*y0 >= 4*self.beadRad*self.beadRad:
                continue

            numTrials += 1
            x = x0 + particleCenter[0]
            y = y0 + particleCenter[1]
            if self.isAvailPoly(particleCenter,(x,y)):
                freePts.append((x,y))
        # bead area = pi r^2 = 1
        # sampling circle area = pi (2r)^2 = 4 pi r^2 = 4
        # free area = len(freePts)/numTrials * sampling circle area
        #           = 4 * len(freePts)/numTrials
        freeArea = 4*len(freePts)/maxTrials
        if freeArea == 0:
            freeArea = 4/maxTrials
        return [particleID,freePts,freeArea]



    def showFreeSpacePoly(self,particleID):
        [pID, freePx] = self.freeSpacePoly(particleID)

        # lots of copied-pasted from showGrid lol
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.scatter(*zip(*self.occupiedPx),marker='.')
        plt.scatter(*zip(*freePx),marker='.',color='red')
        plt.xlim(0,self.gridSize[0])
        plt.ylim(0,self.gridSize[1])
        plt.show()

    def showFreeSpace(self,particleID):
        [pID, freePx] = self.freeSpace(particleID)

        # lots of copied-pasted from showGrid lol
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.scatter(*zip(*self.occupiedPx),marker='.')
        plt.scatter(*zip(*freePx),marker='.',color='red')
        plt.xlim(0,self.gridSize[0])
        plt.ylim(0,self.gridSize[1])
        plt.show()


    def showExclSpace(self,particleCenter):
        exclPx = self.pxExcludedByParticle(particleCenter)

        # lots of copied-pasted from showGrid lol
        fig, ax = plt.subplots()
        for (x,y) in self.particleCenters:
            circ = plt.Circle((x, y), self.beadRad, color='gray', alpha=0.3)
            ax.add_artist(circ)

        plt.scatter(*zip(*self.occupiedPx),marker='.')
        plt.scatter(*zip(*exclPx),marker='.',color='red')

        plt.xlim(0,self.gridSize[0])
        plt.ylim(0,self.gridSize[1])
        plt.show()

    # method options: 1) old px-by-px method; 2) polygons method; 3) MC; 4) analytic
    def entropy(self,makeImg=False,imgFile=None,sbeadFile=None,method=1):
        tic = time.time()

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if sbeadFile == None:
            sbeadFile = nameRoot+'_rad'+str(self.beadRad)+'_Sbead.csv'

        if makeImg:
            # generate an image name, if none provided
            if imgFile == None:
                imgFile = nameRoot+'_rad'+str(self.beadRad)+'_snowflakes.png'

            # initialize a white grid
            scale = 5
            imgArr = np.ones((self.gridSize[1]*scale,self.gridSize[0]*scale,3),dtype=np.uint8) * 255

            # set all the occupied px to blue
            for (x,y) in self.occupiedPx:
                for i in range(scale):
                    for j in range(scale):
                        imgArr[y*scale+j,x*scale+i] = [24,96,148]
            # TODO you should write the img to a csv or smtg
            # TODO you should somehow use showGridNew instead of copy pasting

        S = 0 # total S
        numParts = 0 # number of particles counted

        nbead = len(self.beadShape) # number of px in a particle
        buffer = self.beadRad*2

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            for i in range(len(self.particleCenters)):
                # don't count particles that are too close to the edge
                p = self.particleCenters[i]
                if p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
                p[1] < buffer or p[1] >= self.gridSize[1]-buffer:
                    continue
                
                numParts += 1

                #print("finding freepx for particle at",p)
                #freePx = []
                if method==1:
                    [pID, freePx] = self.freeSpacePoly(i)
                    Si = np.log(len(freePx)/nbead)
                elif method==2:
                    [pID, freePx] = self.freeSpace(i)
                    Si = np.log(len(freePx)/nbead)
                elif method==3:
                    # the second thing returned is rly freePts, a random sampling of the free area
                    # but i am naming it freePx bc i wanna draw it
                    [pID, freePx, freeArea] = self.freeSpaceMC(i)
                    Si = np.log(freeArea)

                S += Si

                if self.usePsi6:
                    psi6 = self.psi6dict[p]
                    writer.writerow([pID,Si,psi6])
                else:
                    writer.writerow([pID,Si])

                if makeImg:
                    # set all the free space px to red
                    if method == 1 or method == 2:
                        for (x,y) in freePx:
                            for i in range(scale):
                                for j in range(scale):
                                    imgArr[y*scale+j,x*scale+i] = [190,25,10]
                    elif method == 3:
                        for (x,y) in freePx:
                            for i in range(scale):
                                for j in range(scale):
                                    imgArr[int(y*scale+j),int(x*scale+i)] = [0,0,0]

        if makeImg:
            # finally time to make & save our image!
            img = PIL.Image.fromarray(imgArr,'RGB')
            img.save(imgFile)

        toc = time.time()

        return [S,numParts,toc-tic]

###########################################################
    def freeSpaceGeoList(self,numProc):
        tic = time.time()

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]

        sbeadFile = nameRoot+'_geo'+'_spaces.csv'


        # pick out the particles to include in entropy calculation
        particlesInGrid = []
        buffer = self.beadRad*2
        for i in range(len(self.particleCenters)):
            # don't count particles that are too close to the edge
            p = self.particleCenters[i]
            if not(p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
            p[1] < buffer or p[1] >= self.gridSize[1]-buffer):
                particlesInGrid.append(i)
        numParts = len(particlesInGrid)

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            # get all snowflakes in parallel
            pool = mp.Pool(numProc)  

            pool_results = []
            pool_results = [pool.apply_async(self.freeSpaceGeo,args=[i]) for i in particlesInGrid]

            pool.close()
            pool.join()

            for r in pool_results:
                [pID,freeA] = r.get()
                
                writer.writerow([pID,freeA])


        toc = time.time()

        return [numParts,toc-tic]
##############################################################################







    def entropyParallel(self,numProc,makeImg=False,imgFile=None,sbeadFile=None,poly=False):
        tic = time.time()

        # generate an Sbead file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if poly:
            polyStr = '_poly'
        else:
            polyStr = '_old' 
            
        if sbeadFile == None:
            sbeadFile = nameRoot+'_rad'+str(self.beadRad)+polyStr+'_Sbead.csv'

        if makeImg:
            # generate an image name, if none provided
            if imgFile == None:
                imgFile = nameRoot+'_rad'+str(self.beadRad)+polyStr+'_snowflakes.png'

            # initialize a white grid
            scale = 5
            imgArr = np.ones((self.gridSize[1]*scale,self.gridSize[0]*scale,3),dtype=np.uint8) * 255

            # set all the occupied px to blue
            for (x,y) in self.occupiedPx:
                for i in range(scale):
                    for j in range(scale):
                        #imgArr[y*scale+j,x*scale+i] = [24,96,148]
                        imgArr[y*scale+j,x*scale+i] = [133,133,133]
            
            cmap = cm.get_cmap('viridis')
            # TODO you should write the img to a csv or smtg
            # TODO you should somehow use showGridNew instead of copy pasting

        S = 0 # total S
        nbead = len(self.beadShape) # number of px in a particle

        # pick out the particles to include in entropy calculation
        particlesInGrid = []
        buffer = self.beadRad*2
        for i in range(len(self.particleCenters)):
            # don't count particles that are too close to the edge
            p = self.particleCenters[i]
            if not(p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
            p[1] < buffer or p[1] >= self.gridSize[1]-buffer):
                particlesInGrid.append(i)
        numParts = len(particlesInGrid)

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            # get all snowflakes in parallel
            pool = mp.Pool(numProc)  

            pool_results = []
            if poly:
                pool_results = [pool.apply_async(self.freeSpacePoly,args=[i]) for i in particlesInGrid]
            else:
                pool_results = [pool.apply_async(self.freeSpace,args=[i]) for i in particlesInGrid]

            pool.close()
            pool.join()

            for r in pool_results:
                [pID,freePx] = r.get()
                Si = np.log(len(freePx)/nbead)
                S += Si
                if self.usePsi6:
                    psi6 = self.psi6dict[p]
                    writer.writerow([pID,Si,psi6])
                else:
                    writer.writerow([pID,Si])

                if makeImg:
                    # set all the free space px to red... or something
                    rgb = [ 255*v for v in cmap(len(freePx)/nbead*3.5)[0:3]]
                    for (x,y) in freePx:
                        for i in range(scale):
                            for j in range(scale):
                                imgArr[y*scale+j,x*scale+i] = rgb
                                #imgArr[y*scale+j,x*scale+i] = [190,25,10]

        if makeImg:
            # finally time to make & save our image!
            img = PIL.Image.fromarray(imgArr,'RGB')
            img.save(imgFile)

        toc = time.time()

        return [S,numParts,toc-tic]


    

    # counts possible configurations of particles in i_neighbors, the lazy way (ie treat each neighbor as if its independent of the others)
    # returns ln(total # configs)
    def neighborConfigs(self,i_neighbors):
        tic = time.time()
        ln_configs = 0

        for i in i_neighbors:
            space = self.freeSpace(i)[1]
            V = len(space)
            ln_configs += np.log(V)
        toc = time.time()
        print('time: '+str(toc-tic))
        return ln_configs

    # TODO i am assuming freePx is sorted
    def numConfigs(self,numBeads,freePx):
        if numBeads == 1:
            return len(freePx)
        
        configs = 0
        n = 0
        ntot = len(freePx) # TODO delete
        for p in freePx:
            if numBeads == 6: # TODO delete
                print(n/ntot) # progress TODO delete
            # consider all the freePx... except the ones that have been already touched by this bead
            newFreePx = freePx[n:]
            # remove all the newly excluded positions from newFreePx
            # TODO i guess when i wrote this i was too edgy to use pxExcludedByParticle. test it properly...
            exclPx = [tuple(map(operator.add, p,x)) for x in self.exclShape]
            #exclPx = pxExcludedByParticle(p)
            for q in exclPx:
                index = bisect.bisect_left(newFreePx, q)
                if index != len(newFreePx) and newFreePx[index] == q:
                    del newFreePx[index]
            configs += self.numConfigs(numBeads-1,newFreePx)
            n += 1

        return configs


    def volumeFraction(self):
        # TODO this would need to change for a hexagonal grid
        return len(self.occupiedPx) / self.gridSize[0] / self.gridSize[1]

    # randomly populates a gridX-by-gridY grid with numBeads beads of radius rad
    # overwrites whatever this PolycrystalGrid used to be
    def populateGridRandomly(self,numBeads,gridX,gridY,rad):

        self.gridSize = [gridX,gridY]
        self.xmin = 0
        self.ymin = 0
        self.beadRad = rad
        self.particleCenters = []
        self.occupiedPx = []
        self.resolution = 1

        self.setBeadShape()
        self.setExclShape()

        i = 0 # number of beads i've placed so far

        tic = time.time()
        while i < numBeads:

            # pick a random position for a new bead
            x = random.randint(self.beadRad,self.gridSize[0]-self.beadRad)
            y = random.randint(self.beadRad,self.gridSize[1]-self.beadRad)

            # if you can, place a bead there
            if self.isAvailable((x,y)):
                self.particleCenters.append((x,y))
                # update occupiedPx
                newPx = self.pxOccupiedByParticle((x,y))
                for p in newPx:
                    bisect.insort_left(self.occupiedPx,p)
                i += 1
            
            
            # if it's been too long, quit
            if time.time()-tic > 60:
                print("timeout: grid production failed. included "+str(i)+" beads")
                break

    def saveGrid(self,fname):
        with open(fname,'w',newline='') as csvFile:
            writer = csv.writer(csvFile,delimiter=',')
            writer.writerow([self.xmin,self.ymin,self.gridSize[0],self.gridSize[1],self.beadRad])
            for (x,y) in self.particleCenters:        
                writer.writerow([x,y])


def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

