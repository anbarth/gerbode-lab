import numpy as np
import csv
import matplotlib.pyplot as plt
import time
import bisect
import random
import operator
import os
import multiprocessing as mp

class PolycrystalGrid:
    # class fields:
    #   xmin: minimum x value
    #   ymin: minimum y value TODO explain these better lol
    #   gridSize: a list that goes [width,height]. note that 0 is included and width and height are excluded.
    #   beadRad: bead radius in px
    #   particleCenters: 
    #   occupiedPx: a list of grid spots occupied by (x,y)s
    #   resolution: 

    # constructor takes the input csv file
    def __init__(self,fname,resolution=1):

        self.resolution = resolution
        self.particleCenters = []
        self.occupiedPx = []
        self.exclShape = []
        self.beadShape = []



        # read in particle centers from csv
        with open(fname) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            first = True
            for row in reader:
                # TODO resolution stuff would probably need to change for a hexagonal grid
                # first line of csv gives xmin, ymin, grid dimensions, and bead radius
                if first:
                    self.xmin = int(row[0])
                    self.ymin = int(row[1])
                    self.gridSize = [round( int(row[2])*resolution ),round( int(row[3])*resolution )]
                    self.beadRad = round( float(row[4])*resolution )
                    first = False
                    continue
                # after the first line, it's particle centers all the way down
                self.particleCenters.append( (round( (float(row[0])-self.xmin)*resolution ),
                                              round( (float(row[1])-self.ymin)*resolution )) )

        # particles always have the same shape -- find that shape
        self.setBeadShape() 
        self.setExclShape() 

        # fill in all the occupied pixels
        for (x0,y0) in self.particleCenters:
            self.occupiedPx.extend(self.pxOccupiedByParticle((x0,y0)))
        self.occupiedPx = sorted(self.occupiedPx)

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
            # TODO tryin out exclusive
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
        pxToCheck = particle[:] # make a COPY!!!!
        for (x,y) in pxToCheck:
            if self.isAvailableIgnore((x,y),particle):
                freePx.append((x,y))
                neighbors = self.getNeighbors((x,y))
                for neighbor in neighbors:
                    if neighbor not in pxToCheck:
                        pxToCheck.append(neighbor)
        return freePx

    def freeSpaceArea(self,particleCenter):
        return len(self.freeSpace(particleCenter))

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

    def entropy(self):
        #print('--- sequential entropy')
        tic = time.time()
        S = 0
        numParts = 0 # number of particles counted
        nbead = len(self.beadShape) # number of px in a particle
        for p in self.particleCenters:
            # don't count particles that are offgrid
            if p[0] < 0 or p[0] >= self.gridSize[0] or p[1] < 0 or p[1] >= self.gridSize[1]:
                continue
            numParts += 1
            nfree = len(self.freeSpace(p)) # number of px available to move to
            S += np.log(nfree/nbead)

        toc = time.time()
        return [S,numParts,toc-tic]
    
    def simpleParallel(self):
        tic = time.time()
        tot = 0
        numParts = len(self.particleCenters)

        with mp.Pool(mp.cpu_count()) as pool:
            pool_results = [pool.apply_async(dist,args=[p,(0,0)]) for p in self.particleCenters]
            #pool_results = pool.imap_unordered(distFrom0, self.particleCenters)
            for r in pool_results:
                dis = r.get()
                tot += dis/numParts

        toc = time.time()
        return [tot,toc-tic]

    def entropyParallel(self,numProc):
        #print('--- parallel entropy')
        tic = time.time()
        
        particlesInGrid = []
        for p in self.particleCenters:
            if p[0] < 0 or p[0] >= self.gridSize[0] or p[1] < 0 or p[1] >= self.gridSize[1]:
                continue
            particlesInGrid.append(p)
        numParts = len(particlesInGrid) # number of particles counted
        #print(numParts)
        nbead = len(self.beadShape) # number of px in a particle
        S = 0

        #with mp.Pool(10) as pool:
        pool = mp.Pool(numProc)
        pool_results = [pool.apply_async(self.freeSpaceArea,args=[p]) for p in particlesInGrid]
        pool.close()
        pool.join()
        #pool_results = pool.imap_unordered(self.freeSpaceArea, particlesInGrid)
        for r in pool_results:
            nfree = r.get() # number of px available to move to
            S += np.log(nfree/nbead)

        toc = time.time()
        return [S,numParts,toc-tic]
        

    def entropyNoEdgeBeads(self):
        S = 0
        numParts = 0
        nbead = len(self.beadShape)
        for p in self.particleCenters:
            # don't count particles that aren't fully on the grid
            if len(self.pxOccupiedByParticle(p)) < nbead:
                continue
            numParts += 1
            nfree = len(self.freeSpace(p)) # number of px available to move to
            S += np.log(nfree/nbead)
        return [S,numParts]

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

