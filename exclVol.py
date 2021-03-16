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
    def __init__(self,fname,rad=0,usePsi6=False):

        self.crystalFile = fname
        self.particleCenters = []
        self.psi6dict = {}
        self.occupiedPx = []
        self.exclShape = []
        self.beadShape = []
        self.snowflakeShape = []
        self.usePsi6 = usePsi6
        self.psi6cutoff = 0.98


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
        # TODO anna you're in the middle of writing this fxn
        # when its done, try running testSimple and get entropy for annasCoolTest10.csv
        for p in self.particleCenters:
            buffer = self.beadRad*2
            if p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
               p[1] < buffer or p[1] >= self.gridSize[1]-buffer:
                continue
            if self.psi6dict[p] >= self.psi6cutoff:
                freePx = self.freeSpace(p)
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

    def freeSpace(self,particleCenter):
        particle = self.pxOccupiedByParticle(particleCenter) # particle includes all (x,y) to ignore
        freePx = []
        #pxToCheck = particle[:] # make a COPY!!!!
        pxToCheck = [particleCenter]
        for (x,y) in pxToCheck:
            if self.isAvailableIgnore((x,y),particle):
                freePx.append((x,y))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)
        return freePx

    def showFreeSpace(self,particleCenter):
        freePx = self.freeSpace(particleCenter)

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

    def entropy(self,makeImg=False,imgFile=None,sbeadFile=None):
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
        Sbead = [] # contribution to S from each bead

        nbead = len(self.beadShape) # number of px in a particle
        buffer = self.beadRad*2

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            for p in self.particleCenters:
                # don't count particles that are too close to the edge
                if p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
                p[1] < buffer or p[1] >= self.gridSize[1]-buffer:
                    continue
                
                numParts += 1
                freePx = self.freeSpace(p)
                Si = np.log(len(freePx)/nbead)
                S += Si

                if self.usePsi6:
                    psi6 = self.psi6dict[p]
                    Sbead.append([Si,psi6])
                    writer.writerow([Si,psi6])
                else:
                    Sbead.append(Si)
                    writer.writerow([Si])

                if makeImg:
                    # set all the free space px to red
                    for (x,y) in freePx:
                        for i in range(scale):
                            for j in range(scale):
                                imgArr[y*scale+j,x*scale+i] = [190,25,10]

        if makeImg:
            # finally time to make & save our image!
            img = PIL.Image.fromarray(imgArr,'RGB')
            img.save(imgFile)

        toc = time.time()

        return [S,numParts,Sbead,toc-tic]

    def entropyParallel(self,numProc,makeImg=False,imgFile=None,sbeadFile=None):
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
                        #imgArr[y*scale+j,x*scale+i] = [24,96,148]
                        imgArr[y*scale+j,x*scale+i] = [133,133,133]
            
            cmap = cm.get_cmap('viridis')
            # TODO you should write the img to a csv or smtg
            # TODO you should somehow use showGridNew instead of copy pasting

        S = 0 # total S
        Sbead = [] # contribution to S from each bead

        nbead = len(self.beadShape) # number of px in a particle

        # pick out the particles to include in entropy calculation
        particlesInGrid = []
        buffer = self.beadRad*2
        for p in self.particleCenters:
            # don't count particles that are too close to the edge
            if not(p[0] < buffer or p[0] >= self.gridSize[0]-buffer or \
            p[1] < buffer or p[1] >= self.gridSize[1]-buffer):
                particlesInGrid.append(p)
        numParts = len(particlesInGrid)

        with open(sbeadFile,'w',newline='') as sbeadFileObj:
            writer = csv.writer(sbeadFileObj)

            nbead = len(self.beadShape) # number of px in a particle
            S = 0
            Sbead = []

            # get all snowflakes in parallel
            pool = mp.Pool(numProc)
            pool_results = [pool.apply(self.freeSpace,args=[p]) for p in particlesInGrid]
            pool.close()
            pool.join()

            for r in pool_results:
                freePx = r.get()
                Si = np.log(len(freePx)/nbead)
                S += Si
                if self.usePsi6:
                    psi6 = self.psi6dict[p]
                    Sbead.append([Si,psi6])
                    writer.writerow([Si,psi6])
                else:
                    Sbead.append(Si)
                    writer.writerow([Si])

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

        return [S,numParts,Sbead,toc-tic]


    

    # counts possible configurations of particles in i_neighbors, the lazy way (ie treat each neighbor as if its independent of the others)
    # returns ln(total # configs)
    def neighborConfigs(self,i_neighbors):
        tic = time.time()
        ln_configs = 0

        for i in i_neighbors:
            p = self.particleCenters[i]
            space = self.freeSpace(p)
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

