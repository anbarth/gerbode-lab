import numpy as np
import csv
import matplotlib.pyplot as plt
import time
import bisect
import random
import operator

#fname = 'myPartsOrderly.csv'
#xmin = 150
#ymin = 350

#fname = 'myParts.csv'
#xmin = 550
#ymin = 350

#beadRad = 8
# grid of pixels includes 0, excludes gridSize
#gridSize = [350,300]

#fname = 'annasNewCrystal.csv'
#fname = 'annasNewCrystalRandom.csv'
#fname = 'small.csv'
#fname = 'empty.csv'
#fname = 'nearestNeighborsHard.csv'
fname = 'testForNearestNeighbors.csv'
i_center = 12
i_neighbors = [6,7,11,13,16,17]
xmin = 0
ymin = 0
beadRad = 0


particleCenters = []
occupiedPx = []

exclShape = []
occShape = []


def populateGridRandomly(numBeads,gridX,gridY):
    global particleCenters
    global occupiedPx

    global gridSize 
    gridSize = [gridX,gridY]

    i = 0
    tic = time.time()
    while i < numBeads:
        x = random.randint(beadRad,gridSize[0]-beadRad)
        y = random.randint(beadRad,gridSize[1]-beadRad)
        if isAvailable((x,y)):
            particleCenters.append((x,y))
            newPx = pxOccupiedByParticle((x,y))
            for p in newPx:
                bisect.insort_left(occupiedPx,p)
            i += 1
        if time.time()-tic > 60:
            print("timeout: grid production failed. included "+str(i)+" beads")
            break

def saveGrid(fname):
    with open(fname,'w',newline='') as csvFile:
        writer = csv.writer(csvFile,delimiter=',')
        writer.writerow([gridSize[0],gridSize[1]])
        for (x,y) in particleCenters:        
            writer.writerow([x,y])


def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

def getNeighbors(p):
    (x,y) = p
    neighbors = []
    if x > 0:
        neighbors.append((x-1,y))
    if y > 0:
        neighbors.append((x,y-1))
    if x < gridSize[0]-1:
        neighbors.append((x+1,y))
    if y < gridSize[1]-1:
        neighbors.append((x,y+1))
    return neighbors

def populateGrid():
    global occupiedPx # necessary bc i have a line of the form "occupiedPx = ..." in this function
    global gridSize

    
    #print('reading in the centers...')
    with open(fname) as csvFile:
        reader = csv.reader(csvFile, delimiter=',')
        first = True
        for row in reader:
            if first:
                gridSize = [int(row[0]),int(row[1])]
                first = False
                continue
            # TODO dont hardcode the offset values here
            particleCenters.append((int(float(row[0])),int(float(row[1])))) #normal
            #particleCenters.append(  (int(float(row[0]))-xmin, gridSize[0]-(int(float(row[1]))-ymin))  ) #myParts
            #particleCenters.append(  (int(float(row[0]))-150, gridSize[0]-(int(float(row[1]))-350))  ) #myPartsOrderly

    setOccShape()
    setExclShape()

    #print('filling in occupied px...')
    for (x0,y0) in particleCenters:
        occupiedPx.extend(pxOccupiedByParticle((x0,y0)))
    occupiedPx = sorted(occupiedPx)


    
    
def showGrid():
    fig, ax = plt.subplots()
    for (x,y) in particleCenters:
        circ = plt.Circle((x, y), beadRad, color='gray', alpha=0.3)
        ax.add_artist(circ)

    '''circ = plt.Circle(particleCenters[i_center], beadRad, color='purple', alpha=0.5)
    ax.add_artist(circ)

    for x in i_neighbors:
        circ = plt.Circle(particleCenters[x], beadRad, color='green', alpha=0.5)
        ax.add_artist(circ)'''

    plt.scatter(*zip(*occupiedPx),marker='.')
    plt.xlim(0,gridSize[0])
    plt.ylim(0,gridSize[1])
    #plt.xticks(np.arange(0, gridSize[0], step=1))
    #plt.yticks(np.arange(0, gridSize[1], step=1))
    #plt.grid(b=True,which='both',axis='both')
    plt.show()

# returns the shape of a generic excluded volume
# requires a seed position x0,y0 that's far from the edges
# i've chosen the center of the grid
# returns as a list
# to get the excluded volume for a bead centered at (x,y), add (x,y) to every element of the list
# TODO is it a problem that that ^ gives you positions that aren't actually on the grid??
def setExclShape():
    global exclShape

    x0 = int(gridSize[0]/2)
    y0 = int(gridSize[1]/2)
    # px excluded by p are all the positions that conflict with a bead at p
    # get all the spots occupied by a bead at p
    # sort that shit bc we're gonna be searching it a LOT
    pxList = sorted(pxOccupiedByParticle((x0,y0)))

    # get ready to check for excluded positions, starting at the center p
    pxToCheck = [(x0,y0)]
    excluded = []

    for (x,y) in pxToCheck:
        if conflict((x,y),pxList):
            excluded.append((x-x0,y-y0))
            neighbors = getNeighbors((x,y))
            for neighbor in neighbors:
                if neighbor not in pxToCheck:
                    pxToCheck.append(neighbor)

    exclShape = excluded
    return excluded


def setOccShape():
    global occShape

    x0 = int(gridSize[0]/2)
    y0 = int(gridSize[1]/2)

    pxToCheck = [(x0,y0)]
    occupied = []

    for (x,y) in pxToCheck:
        if dist((x0,y0),(x,y)) <= beadRad:
            occupied.append((x-x0,y-y0))
            neighbors = getNeighbors((x,y))
            for neighbor in neighbors:
                if neighbor not in pxToCheck:
                    pxToCheck.append(neighbor)

    occShape = occupied
    return occupied


def pxExcludedByParticle(p):

    # px excluded by p are all the positions that conflict with a bead at p

    # get all the spots occupied by a bead at p
    # sort that shit bc we're gonna be searching it a LOT
    pxList = sorted(pxOccupiedByParticle(p))

    # get ready to check for excluded positions, starting at the center p
    pxToCheck = [p]
    excluded = []

    for (x,y) in pxToCheck:
        if conflict((x,y),pxList):
            excluded.append((x,y))
            neighbors = getNeighbors((x,y))
            for neighbor in neighbors:
                if neighbor not in pxToCheck:
                    pxToCheck.append(neighbor)


    return excluded

def pxOccupiedByParticleNew(p):
    #return [tuple(map(operator.add, p,x)) for x in occShape]
    pxOcc = []
    for (x,y) in occShape:
        xpos = x+p[0]
        ypos = y+p[1]
        if xpos >= 0 and xpos < gridSize[0] and ypos >= 0 and ypos < gridSize[1]:
            pxOcc.append((xpos,ypos))
    return pxOcc


def pxOccupiedByParticle(p):
    (x0,y0) = p
    pxToCheck = [(x0,y0)]
    occupied = []

    for (x,y) in pxToCheck:
        if dist((x0,y0),(x,y)) <= beadRad:
            occupied.append((x,y))
            neighbors = getNeighbors((x,y))
            for neighbor in neighbors:
                if neighbor not in pxToCheck:
                    pxToCheck.append(neighbor)

    return occupied

# TODO isAvailable is just this but with pxList=occupiedPx
# p=(x0,y0) is a test position. does it conflict with the positions in pxList?
def conflict(p,pxList):
    (x0,y0) = p
    # pretend to put a bead at (x0,y0)
    dangerZone = pxOccupiedByParticle((x0,y0))
    # now see if any of the positions occupied by the bead are in pxList
    for (x,y) in dangerZone:
        index = bisect.bisect_left(pxList,(x,y))
        if index != len(pxList) and pxList[index] == (x,y):
        #if (x,y) in pxList:
            return True
    return False

def isAvailable(p):
    (x0,y0) = p
    dangerZone = pxOccupiedByParticle((x0,y0))
    for (x,y) in dangerZone:
        index = bisect.bisect_left(occupiedPx,(x,y))
        if index != len(occupiedPx) and occupiedPx[index] == (x,y):
        #if (x,y) in occupiedPx:
            return False
    return True

def isAvailableIgnore(p,particle):
    (x0,y0) = p
    dangerZone = pxOccupiedByParticle((x0,y0))
    for (x,y) in dangerZone:
        #if (x,y) in occupiedPx and (x,y) not in particle:
        if (x,y) not in particle:
            index = bisect.bisect_left(occupiedPx,(x,y))
            if index != len(occupiedPx) and occupiedPx[index] == (x,y):
                return False
    return True

def freeSpace(particleCenter):
    tic = time.time()
    particle = pxOccupiedByParticle(particleCenter) # particle includes all (x,y) to ignore
    freePx = []
    pxToCheck = particle[:] # make a COPY!!!!
    for (x,y) in pxToCheck:
        if isAvailableIgnore((x,y),particle):
            freePx.append((x,y))
            neighbors = getNeighbors((x,y))
            for neighbor in neighbors:
                if neighbor not in pxToCheck:
                    pxToCheck.append(neighbor)

    
    toc = time.time()
    #print("time: "+str(toc-tic))
    #print("px checked: "+str(len(pxToCheck)))
    return freePx

def showFreeSpace(particleCenter):
    freePx = freeSpace(particleCenter)

    # lots of copied-pasted from showGrid lol
    fig, ax = plt.subplots()
    for (x,y) in particleCenters:
        circ = plt.Circle((x, y), beadRad, color='gray', alpha=0.3)
        ax.add_artist(circ)

    plt.scatter(*zip(*occupiedPx),marker='.')
    plt.scatter(*zip(*freePx),marker='.',color='red')

    plt.xlim(0,gridSize[0])
    plt.ylim(0,gridSize[1])
    #plt.xticks(np.arange(0, gridSize[0], step=1))
    #plt.yticks(np.arange(0, gridSize[1], step=1))
    #plt.grid(b=True,which='both',axis='both')

    plt.show()


def showExclSpace(particleCenter):
    #shadow = exclShape()
    exclPx = [tuple(map(operator.add, particleCenter,x)) for x in exclShape]

    # lots of copied-pasted from showGrid lol
    fig, ax = plt.subplots()
    for (x,y) in particleCenters:
        circ = plt.Circle((x, y), beadRad, color='gray', alpha=0.3)
        ax.add_artist(circ)

    plt.scatter(*zip(*occupiedPx),marker='.')
    plt.scatter(*zip(*exclPx),marker='.',color='red')

    plt.xlim(0,gridSize[0])
    plt.ylim(0,gridSize[1])
    #plt.xticks(np.arange(0, gridSize[0], step=1))
    #plt.yticks(np.arange(0, gridSize[1], step=1))
    #plt.grid(b=True,which='both',axis='both')

    plt.show()

def entropy():
    tic = time.time()
    S = 0
    n = 0
    for p in particleCenters:
        #print(n)
        n += 1
        space = freeSpace(p)
        V = len(space)
        particleVol = len(pxOccupiedByParticle(p))
        S += np.log(V/particleVol)
        #S += np.log(V)
    toc = time.time()
    print('time: '+str(toc-tic))
    return S

# counts possible configurations of particles in i_neighbors, the lazy way (ie treat each neighbor as if its independent of the others)
# returns ln(total # configs)
def neighborConfigs():
    tic = time.time()
    ln_configs = 0

    for i in i_neighbors:
        p = particleCenters[i]
        space = freeSpace(p)
        V = len(space)
        ln_configs += np.log(V)
    toc = time.time()
    print('time: '+str(toc-tic))
    return ln_configs

# TODO i am assuming freePx is sorted
def numConfigs(numBeads,freePx):
    if numBeads == 1:
        return len(freePx)
    
    configs = 0
    n = 0
    for p in freePx:
        #print(n)
        # consider all the freePx... except the ones that have been already touched by this bead
        newFreePx = freePx[n:]
        # remove all the newly excluded positions from newFreePx
        exclPx = [tuple(map(operator.add, p,x)) for x in exclShape]
        #exclPx = pxExcludedByParticle(p)
        for q in exclPx:
            index = bisect.bisect_left(newFreePx, q)
            if index != len(newFreePx) and newFreePx[index] == q:
                del newFreePx[index]

        configs += numConfigs(numBeads-1,newFreePx)
        n += 1

    return configs


def volumeFraction():
    return len(occupiedPx) / gridSize[0] / gridSize[1]

populateGrid()
#populateGridRandomly(36)
#print(len(particleCenters))
#print(volumeFraction())
#print(len(pxOccupiedByParticle(particleCenters[152])))
#print(neighborConfigs())

showGrid()
#print(pxExcludedByParticle((2,2)))
showFreeSpace((2,2))
#showFreeSpace(particleCenters[5])
#pxList = pxOccupiedByParticle((21,24))
#print(conflict((16,24),pxList))
#saveGrid('test.csv')
#showFreeSpace((21,14))
#showFreeSpace((21,32))

#totFreeSpace = sorted( freeSpace((21,14)) + freeSpace((21,32)) )
#tic = time.time()
#print(numConfigs(3,totFreeSpace))
#print('time: '+str(time.time()-tic))



'''i = 100
while i < 103:
    print(i)
    print(particleCenters[i])
    showFreeSpace(particleCenters[i])
    i += 1'''