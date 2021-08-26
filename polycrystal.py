import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import multiprocessing as mp
import myGeo
import importlib
import os
from matplotlib import path
importlib.reload(myGeo)
import random

####   welcome to the polycrystal class!
####   the primary purpose of this class is to find the entropy associated with a polycrystal
####   by finding the free space available to each particle.
####   more details on how to use this class can be found on the gerbode lab website.
####   anna barth 2021

class Polycrystal:
    # class fields:
    #   crystalFile: the name of the csv containing a list of particle positions
    #   beadRad: the particle radius
    #   particleCenters: a list of particle positions, stored as tuples,
    #        indexed by particle ID-1 (i.e. particle 1's position is particleCenters[0]).
    #   neighbs: a dictionary of each particle's nearest neighbors, stored as a list of 
    #        particle IDs and sorted in CCW order, indexed by particle ID
    #        (i.e. particle 1's nearest neighbors are neighbs[1]).
    #   windowVertices: the vertices that define the window. particles whose centers are inside this polygon
    #        will be counted when calculating entropy (unless windowOverride is on).
    #        stored as a list of (x,y) tuples in CW or CCW order.
    #   windowOverride: a boolean OR a string. if True, then we will decide which particles get counted
    #        by reading off the 4th column of the crystal csv (rather than looking at which
    #        particles fall within a certain window). if a string, it's interpreted as the name
    #        of a csv file that says, line-by-line, which particles should get counted.
    #   countParticle: a list of bools, indexed by particle ID-1, that says whether each particle
    #        should be counted when calculating entropy
    #   displayWindow: the limits of the window that should actually be displayed when making
    #        a picture of the polycrystal. stored as a list [left, right, bot, top], and determined based
    #        on windowVertices.

    # ===================================================================================================   
    # =========================================== CONSTRUCTOR ===========================================
    # ===================================================================================================   

    # input fname: the crystal csv file name
    # optional input neighbFile: the nearest neighbor file name
    #      (if not given, i will assume it's at crysNeighbs/crystalFileName_neighbs.csv)
    # optinal input windowOverride: see description above
    # optinal input radius: if provided, i'll use this instead of the value in the crystal csv
    def __init__(self,fname,neighbFile=None,windowOverride=False,radius=0):

        self.crystalFile = fname
        self.particleCenters = []
        self.neighbs = {}
        self.windowVertices = []
        self.displayWindow = [10000,-10000,10000,-1000] # dummy values to start
        self.countParticle = []

        # read the crystal csv file
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
                # second line gives polygon window as [x1,y1,x2,y2...]
                if lineNum == 2:                  
                    if len(row) % 2 != 0:
                        print("warning: row 2 should have an even number of entries, it doesn't")
                    for i in range(int(len(row)/2)):
                        myX = float(row[2*i])
                        myY = float(row[2*i+1])
                        self.windowVertices.append((myX,myY))

                        # adjust display window bounds:
                        # left
                        if myX < self.displayWindow[0]:
                            self.displayWindow[0] = myX
                        # right
                        if myX > self.displayWindow[1]:
                            self.displayWindow[1] = myX
                        # bot
                        if myY < self.displayWindow[2]:
                            self.displayWindow[2] = myY
                        # top
                        if myY > self.displayWindow[3]:
                            self.displayWindow[3] = myY
                    lineNum += 1
                    continue
                        

                # after the first line, it's particle centers all the way down
                # col 0 gives particle ID, but that's just so that the csv is human-readable.
                # it's redundant information bc the particles need to be in particle ID order anyway

                # col 1 and 2 are (x,y) position
                p = ( (float(row[1])) , 
                      (float(row[2])) )
                self.particleCenters.append(p)
                
                # if windowOverride is on, col 3 says if this particle should be counted when calculating S
                if windowOverride == True:
                    self.countParticle.append(bool(int(row[3])))

        # if windowOverride is off, decide which particles to count based on if they're inside the window
        if windowOverride == False:
            windowPath = path.Path(self.windowVertices)
            self.countParticle = windowPath.contains_points(self.particleCenters)

        # if windowOverride is a string rather than a bool, it's a csv file that specifies who's in and who's out
        if isinstance(windowOverride,str):
            with open(windowOverride) as windowFile:
                wreader = csv.reader(windowFile, delimiter=',')
                for wrow in wreader:
                    self.countParticle.append(bool(int(wrow[1])))


        # get ready to read the neighbors file
        # if no neighbor file is given, then look for the neighbor file with corresponding name
        # e.g. if the crystal file is crystals/greg.csv, look for crysNeighbs/greg_neighbs.csv
        if neighbFile == None:
            # chop off any part of the name before the first slash
            # (bc crystalFile usually starts with "crystals/")
            i = self.crystalFile.find('/') 
            nameRoot = self.crystalFile[i+1:-4]
            neighbFile = "crysNeighbs/"+nameRoot+"_neighbs.csv"
        else:
            neighbFile = "crysNeighbs/"+neighbFile

        # time to read the neighbors in!
        with open(neighbFile) as csvFile:
            reader = csv.reader(csvFile, delimiter=',')
            for row in reader:
   
                # column 0 gives the particles whomst neighbors we will see
                partID = int(row[0])
                
                # subsequent elements are that particle's neighbors
                myNeighbs = []
                for i in range(1,len(row)):
                    myNeighbs.append(int(row[i]))


                # sortKey is a function that returns a point's angle and distance relative to this particle
                sortKey = myGeo.make_clockwiseangle_and_distance(self.particleCenters[partID-1])

                # use sortKey to sort the neighbors in cw order, which is convenient in freeSpace
                sortMe = np.array([sortKey(self.particleCenters[nnID-1]) for nnID in myNeighbs],dtype="f,f")
                sortedIndices = np.argsort(sortMe)
                self.neighbs[partID] = [myNeighbs[i] for i in sortedIndices]
           

    # ===================================================================================================            
    # ======================================== DISPLAY FUNCTIONS ========================================
    # ===================================================================================================   

    # display the polycrystal
    def show(self):
        self.showHighlight([])

    # display the polycrystal and highlight some particles (specified by pIDs) in red
    def showHighlight(self,pIDs):

        fig, ax = plt.subplots()
        ax.set_aspect(1)

        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]

            if i+1 in pIDs:
                myFaceColor = '#5c0000'
            else:
                myFaceColor = 'gray'

            # draw the particle
            circ = plt.Circle(p, self.beadRad, facecolor=myFaceColor,edgecolor=None, alpha=0.5)
            ax.add_artist(circ)
    
            # if this particle is in the window, add a small dot on its center
            if self.countParticle[i]:
                circ = plt.Circle(p, self.beadRad/5, facecolor='k', edgecolor=None)
                ax.add_artist(circ)
                # sometimes it's also nice to display pIDs
                #plt.text(p[0]+self.beadRad/15,p[1]+self.beadRad/15,str(i+1))

        plt.xlim(self.displayWindow[0]-2*self.beadRad,self.displayWindow[1]+2*self.beadRad)
        plt.ylim(self.displayWindow[2]-2*self.beadRad,self.displayWindow[3]+2*self.beadRad)

        plt.show()
        # you can also automatically save the image if you want
        #strPos = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        #nameRoot = self.crystalFile[strPos+1:-4]
        #fig.savefig(nameRoot+"_img.png",dpi=900)
        plt.close(fig)

    # a silly function to display grainsplitting sims, just to help figure-making
    # idk you can like, delete this function
    def showGrainSplitSim(self,R_grain,isInitialFrame):
        x_max = 2*R_grain+200

        fig, ax = plt.subplots()
        ax.set_aspect(1)

        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]

            if self.countParticle[i] and isInitialFrame:
                myFill = True
                myOrder = 5
                faceColor = '#FFFFFF'
                edgeColor = '#6b6b6b'
            elif p[0] <= x_max/2: #or i==993 or i==954 or i==995 or i==919:
                myFill = True
                myOrder = 0
                faceColor = '#bababa'
                edgeColor = '#bababa'
            else:
                myFill = True
                myOrder = 0
                faceColor = '#636363'
                edgeColor = '#636363'

            circ = plt.Circle(p, self.beadRad, facecolor=faceColor,edgecolor=edgeColor, alpha=1,fill=myFill,zorder=myOrder)
            ax.add_artist(circ)

        plt.xlim(self.displayWindow[0]-4*self.beadRad,self.displayWindow[1]+4*self.beadRad)
        plt.ylim(self.displayWindow[2]-4*self.beadRad,self.displayWindow[3]+4*self.beadRad)
        plt.show()
        plt.close(fig)

    # display a zoomed-in picture of a particle with red circles of radius 2R
    # around its neighbors. very useful for diagnosing issues with freeSpace
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
            #plt.text(q[0]+self.beadRad/15,q[1]+self.beadRad/15,str(i+1))

        p = self.particleCenters[pID-1]
        #plt.text(p[0]+self.beadRad/15,p[1]+self.beadRad/15,str(pID))

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

    # display the free space of the particle specified by particleID
    def showFreeSpace(self,particleID):
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

        plt.plot(freeSpaceCurveX,freeSpaceCurveY,'r')
        plt.xlim(min(freeSpaceCurveX)-self.beadRad/3,max(freeSpaceCurveX)+self.beadRad/3)
        plt.ylim(min(freeSpaceCurveY)-self.beadRad/3,max(freeSpaceCurveY)+self.beadRad/3)
        plt.show()
        plt.close(fig)
        return freeArea

    # ===================================================================================================            
    # =========================================== FREE SPACE ============================================
    # ===================================================================================================   
    
    # find the shape & area of free space associated with the particle specified by particleID
    # returns a list containing: 0. particleID (returning this makes it possible to parallelize entropy, 
    #                               which i don't do anymore lol)
    #                            1. area of the free space (normalized)
    #                            2. array of x-values specifying the shape of the free space
    #                            3. array of y-values specifying the shape of the free space
    # optional input imgMaking: if True, i will automatically save an image of this particle's free
    #                             space if/when something weird happens. the image will be saved
    #                             to badFreeSpace/path/to/crystal/crystalFileName_freespace_pID.png
    def freeSpace(self,particleID,imgMaking=False):
        # a bool to keep track of whether anything weird has happened
        thisWentSmoothly = True

        # return 0 for particles that shouldn't be counted
        if not self.countParticle[particleID-1]:
            return 0

        center = self.particleCenters[particleID-1]
        nnIDs = self.neighbs[particleID]
        nns = [self.particleCenters[nnID-1] for nnID in nnIDs]

        # set up lists to keep track of crossing points:
        # (x,y) for crossing points between neighbors' excl area circles
        crossingPts = [] 
        # (circle1,circle2) specifying, for each entry in crossingPts, which neighbors' circles are crossing
        # stored as the neighbor's index in the list nnIDs
        crossingPairs = [] 

        # go over all pairs of neighbors
        for i in range(len(nns)):
            for j in range(i+1,len(nns)): 
                # find crossing points
                myCrossingPts = myGeo.circIntersections(nns[i][0], nns[i][1], 2*self.beadRad, \
                                                        nns[j][0], nns[j][1], 2*self.beadRad)
                
                # no crossing points? skip
                if myCrossingPts == None:
                    continue
                
                # keep only the closest crossing point
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
                
                # if you've made it this far, add to our list of crossing points
                if keepMe:
                    crossingPts.append(closestCrossingPt)
                    crossingPairs.append( (i,j) )
        
        # check for nearest neighbor network issue #1: neighbors that have >2 crossing points
        # fix it by keeping only the closest 2
        for i in range(len(nns)):
            # find the crossing points on this circle
            myPtInds = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPtInds.append(j)

            # if there's <=2 crossing points, there's no issue
            if len(myPtInds) <= 2:
                continue
            
            if thisWentSmoothly:
                print(particleID,"has extraneous neighbors")
                thisWentSmoothly = False

            # find the two crossing points closest to the central particle
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
                    crossingPts = np.delete(crossingPts,j,0)
                    crossingPairs = np.delete(crossingPairs,j,0)

        # check for nearest neighbor network issue #2: neighbors that have only 1 crossing point
        # fix it by deleting the straggler crossing point
        for i in range(len(nns)):
            # find the crossing points on this circle
            myPtInds = []
            for j in range(len(crossingPts)):
                if crossingPairs[j][0] == i or crossingPairs[j][1] == i:
                    myPtInds.append(j)

            if len(myPtInds) == 1:
                j = myPtInds[0]
                crossingPts = np.delete(crossingPts,j,0)
                crossingPairs = np.delete(crossingPairs,j,0)

        # oh no! something has gone horribly wrong. bail out
        if len(crossingPts) == 0:
            print("free area appears to be 0 or negative for "+str(particleID)+", go check it out")
            return [particleID,0,[center[0]],[center[1]]]

        # find the area of the polygon bounded by the crossing points
        # first, find a point (x_inside, y_inside) that is inside the polygon
        [x_inside,y_inside] = myGeo.centroid([crossPt[0] for crossPt in crossingPts],[crossPt[1] for crossPt in crossingPts])
        # now define sortKey, a function that returns a point's angle and distance relative to (x_inside,y_inside)
        # this allows us to sort the crossing points in ccw order, which is a necessary pre-req for polyArea
        sortKey = myGeo.make_clockwiseangle_and_distance((x_inside,y_inside))
        myArea = myGeo.polyArea(sorted(crossingPts,key=sortKey))

        # find the segments of excluded area that protrude into the polygon
        # these lists will store the (x,y) points that define the free space's perimeter
        freeSpaceCurveX = []
        freeSpaceCurveY = []

        # iterate over the excluded area circles and cut out the appropriate segment for each
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
                print(particleID,"somehow still has an issue, come fix it")
                thisWentSmoothly = False

            # vectors pointing from this neighbor to each crossing pt
            vec1 = (myPts[0][0]-nns[i][0],myPts[0][1]-nns[i][1])
            vec2 = (myPts[1][0]-nns[i][0],myPts[1][1]-nns[i][1])

            # angle of each vector, relative to the positive x axis
            theta1 = np.arctan2(vec1[1],vec1[0])
            theta2 = np.arctan2(vec2[1],vec2[0])

            # arctan2's range is [-pi,pi], meaning there's an awkward jump at pi.
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
            # segment area = (1/2) * (theta-sin(theta)) * (2R)^2
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
            print(particleID,"has weirdly large free area, worth checking out")
            thisWentSmoothly = False

        # if anything fishy happened and imgMaking is on, make & save a picture of the free space
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

            plt.plot(freeSpaceCurveX,freeSpaceCurveY,'r')
            plt.xlim(min(freeSpaceCurveX)-self.beadRad/3,max(freeSpaceCurveX)+self.beadRad/3)
            plt.ylim(min(freeSpaceCurveY)-self.beadRad/3,max(freeSpaceCurveY)+self.beadRad/3)

            # generate a name for the image file
            # chop off the "crystals/" prefix from the crystal file name
            strPos = self.crystalFile.find('/') 
            nameRoot = self.crystalFile[strPos+1:-4]
            myFileName = "badFreeSpace/"+nameRoot+"_freespace_"+str(particleID)+".png"
            try:
                fig.savefig(myFileName, dpi=900)
            # if the folder doesn't exist yet, make it
            except FileNotFoundError:
                strPos = nameRoot.rfind('/')
                os.makedirs("badFreeSpace/"+nameRoot[0:strPos])
                fig.savefig(myFileName, dpi=900)
            plt.close(fig)

        # collapse onto the finish line of this behemoth function
        return [particleID,myArea/(np.pi*self.beadRad*self.beadRad),freeSpaceCurveX,freeSpaceCurveY]


    # ===================================================================================================            
    # ============================================= ENTROPY =============================================
    # ===================================================================================================   

    # find the dimensionless entropy S = \sum_i \ln\frac{v_i}{\pi R^2}
    # returns S, and creates 2 files:
    # (1) a csv where each row is [particle ID, free space]. by default saved to freeSpaceCSVs/crystalFileName_freeSpaces.csv
    # (2) an image showing each particle's free space. by default saved to freeSpaceDiagrams/crystalFileName_snowflakes.png
    # optional inputs freeSpaceFile and imgFile allow you to specify
    # the names of these files rather than using the default
    def entropy(self,freeSpaceFile=None,imgFile=None):
        # generate file names, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if freeSpaceFile == None:
            freeSpaceFile = 'freeSpaceCSVs/'+nameRoot+'_freeSpaces.csv'
        if imgFile == None:
            imgFile = 'freeSpaceDiagrams/'+nameRoot+'_snowflakes.png'

        # setting up the image
        fig, ax = plt.subplots()
        #ax.set_facecolor([0.25,0.25,0.25])
        cmap = cm.get_cmap('viridis')
        plt.xlim(self.displayWindow[0]-self.beadRad,self.displayWindow[1]+self.beadRad)
        plt.ylim(self.displayWindow[2]-self.beadRad,self.displayWindow[3]+self.beadRad)
        ax.set_aspect(1) # this makes it square
        plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False) # turn off the ticks
        plt.tick_params(axis='y',which='both',left=False,right=False,labelleft=False) # turn off the ticks

        S = 0 # total S

        with open(freeSpaceFile,'w',newline='') as freeSpaceFileObj:
            writer = csv.writer(freeSpaceFileObj)

            # iterate over the particles
            for i in range(len(self.particleCenters)):
                p = self.particleCenters[i]

                # draw this particle
                circ = plt.Circle(p, self.beadRad, facecolor=[0.6,0.6,0.6],edgecolor='black',linewidth=0,alpha=1,zorder=0)
                #circ = plt.Circle(p, self.beadRad, facecolor='#b8b8b8',edgecolor='black',linewidth=0,alpha=1,zorder=0)
                ax.add_artist(circ)

                # if this particle shouldn't be counted, then stop now
                if not self.countParticle[i]:
                    continue

                # draw a black dot to indicate this particle is included
                circ = plt.Circle(p, self.beadRad/15, facecolor='k', edgecolor=None,zorder=10)
                ax.add_artist(circ) 

                # find the free area -- note that particleID = i+1 !!
                (pID,freeArea,freeSpaceCurveX,freeSpaceCurveY) = self.freeSpace(i+1,imgMaking=True)

                if freeArea <= 0:
                    # if freeArea<=0, something went very wrong!
                    # this keeps it from crashing, but you should really
                    # go back and fix the issue
                    Si = 0
                else:
                    Si = np.log(freeArea)
                S += Si

                # write down this particle's free area
                writer.writerow([pID,freeArea])

                # draw this particle's free area

                # cmap takes a number in [0,1) to a color
                rgb = cmap((freeArea)*25)[0:3]
                #rgb = cmap((freeArea-0.035)*6.9)[0:3]
                #ax.fill(freeSpaceCurveX,freeSpaceCurveY, facecolor=rgb,edgecolor='black',lw=0.15)
                #ax.fill(freeSpaceCurveX,freeSpaceCurveY, facecolor=[0.8,0.8,0.8],edgecolor='black',lw=0.15)

                # use this code to plot free space at 2x its size
                cen = myGeo.centroid(freeSpaceCurveX,freeSpaceCurveY)
                biggerCurveX = 2*(np.array(freeSpaceCurveX)-cen[0])+cen[0]
                biggerCurveY = 2*(np.array(freeSpaceCurveY)-cen[1])+cen[1]
                ax.fill(biggerCurveX,biggerCurveY, facecolor=rgb,edgecolor='black',lw=0.15,zorder=5)
        

        fig.savefig(imgFile, dpi=900)#, transparent=True)
        plt.close(fig)
        return S

    # an attempt to find entropy by using a monte carlo method to 
    # directly measure \Omega, the volume of accessible, allowed configurations
    # i never got \Omega to be nonzero tho. i think you'd need to run a lot of trials.
    def entropyMC(self,numTrials):
        successfulTrials = 0
        for i in range(numTrials):
            successfulTrials += self.MCtrial()
        return successfulTrials/numTrials

    # run one trial for entropyMC
    def MCtrial(self):
        # randomly displace each particle, put the new positions in trialCenters
        trialCenters = []
        for i in range(len(self.particleCenters)):
            p = self.particleCenters[i]

            # x0 and y0 are in [-2beadRad,+2beadRad)
            x0 = random.random()*4*self.beadRad - 2*self.beadRad
            y0 = random.random()*4*self.beadRad - 2*self.beadRad
            
            trialCenters.append( (p[0]+x0, p[1]+y0) )
        
        # look for overlap between the new particle positions
        for pID in range(1,len(trialCenters)+1):
            # look for overlap only with nearest neighbors
            myNeighbs = self.neighbs[pID]
            for nnID in myNeighbs:
                # if there are overlapping particles, don't count this trial
                if dist(trialCenters[pID-1],trialCenters[nnID-1]) < self.beadRad*2:
                    return 0
        
        # no overlapping particles! count this trial!
        return 1

    # find the dimensionless entropy S by computing all particles' free areas in parallel
    # this function is outdated, i haven't used it in a long time and i don't even know if it would work anymore...
    # returns S and creates a csv where each row is [particle ID, free space],
    # by default saved to freeSpaceCSVs/crystalFileName_freeSpaces.csv
    # optional input freeSpaceFile allows you to specify the name of this file rather than using the default
    # input numProc specifies the number of processes for the parallel stuff
    def entropyParallel(self,numProc,freeSpaceFile=None):

        # generate free space file name, if none provided
        i = self.crystalFile.rfind('/') # chop off any part of the name before a slash
        nameRoot = self.crystalFile[i+1:-4]
        if freeSpaceFile == None:
            freeSpaceFile = 'freeSpaceCSVs/'+nameRoot+'_freeSpaces.csv'

        S = 0 # total S

        # pick out the particles to include in entropy calculation
        particlesInGrid = []
        for i in range(len(self.particleCenters)):
            # don't count particles that are outside the window
            if self.countParticle[i]:
                particlesInGrid.append(i+1)

        with open(freeSpaceFile,'w',newline='') as freeSpaceFileObj:
            writer = csv.writer(freeSpaceFileObj)

            # get all snowflakes in parallel
            pool = mp.Pool(numProc) 

            # apply_async does not necessarily return results in the order you gave the inputs
            # that's why i made entropy return particleID, so that you can still tell
            # which particle goes with which free area
            pool_results = [pool.apply_async(self.freeSpace,args=[pID]) for pID in particlesInGrid]
            pool.close()
            pool.join()

            for r in pool_results:
                [pID,freeArea,freeSpaceX,freeSpaceY] = r.get()
                Si = np.log(freeArea)
                S += Si
                writer.writerow([pID,Si])

        return S

    # ===================================================================================================            
    # ============================================== MISC ===============================================
    # ===================================================================================================  

    # returns the approximate volume fraction within the window
    # assumes area occupied = N pi r^2
    # this isn't accounting for particles that fall partially outside the window
    # i wrote a function that does that, but it's in matlab RIP
    def areaFraction(self):  
        windowPath = path.Path(self.windowVertices)
        inWindow = windowPath.contains_points(self.particleCenters)
        numParts = sum(inWindow)
        areaOccupied = numParts*np.pi*self.beadRad*self.beadRad
        totalArea = myGeo.polyArea(self.windowVertices)
        return areaOccupied/totalArea
    
    # returns the smallest distance between any two particles in the whole polycrystal
    def minDist(self):
        minD = self.beadRad*10 # dummy value to start
        minPair = (-1,-1)
        for i in range(len(self.particleCenters)):
            pID = i+1
            nnIDs = self.neighbs[pID]
            for nnID in nnIDs:
                myD = dist(self.particleCenters[pID-1],self.particleCenters[nnID-1])
                if myD <= minD:
                    minD = myD
                    minPair = (pID,nnID)
        return (minD,minPair)

    # returns the number of particles that get counted when computing entropy
    def numParts(self):
        return sum(self.countParticle)


def dist(p1,p2):
    (x1,y1) = p1
    (x2,y2) = p2
    return np.sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2))

