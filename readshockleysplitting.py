import numpy as np
import matplotlib.pyplot as plt
import csv


# ancient attempt to estimate the energetic cost of grain splitting
# by using read shockely. this does not remotely work lol


class Crystal(object):

    def __init__(self, granule=[], phi_f=[], boundary=[]):
        ''' granule should be a list of each granule's initial orientation angle
            phi_f should be a list of each granule's final orientation angle
            boundary should be a list of GBs, expressed as [length, index of grain 1, index of grain 2]
        '''
        if len(granule) != len(phi_f):
            raise Exception('granules and phi_f must be same length')
        # copy element-by-element to avoid fucking around w external variables
        self._granule = [x for x in granule]
        self._phi_f = [x for x in phi_f]
        self._boundary = [x for x in boundary]
        self._deltaPhi = [ phi_f[i] - granule[i] for i in range(len(granule)) ]
    
    def totalEnergy(self):
        energy = 0
        for GB in self._boundary:
            '''theta = np.abs( self._granule[GB[1]] - self._granule[GB[2]] )
            if theta == 0:
                continue
            energy += theta * (1 - np.log(theta)) * GB[0]'''
            energy += self.gbEnergy(GB)
        return energy

    def granuleEnergy(self,granNum):
        energy = 0
        for GB in self._boundary:
            if GB[1] == granNum or GB[2] == granNum:
                energy += self.gbEnergy(GB)
        return energy

    def gbEnergyByNum(self,GBnum):
        return self.gbEnergy(self._boundary[GBnum])
        
    def gbEnergy(self,GB):
        theta = np.abs( self._granule[GB[1]] - self._granule[GB[2]] )
        if theta == 0:
            return 0
        return theta * (1 - np.log(theta)) * GB[0]

    


    def advanceInTime(self,dt):
        # rotate all granules
        for i in range(len(self._granule)):
            #if self._granule[i] >= self._phi_f[i]:
                # already reached the final time
            #    break
            # total time is always 100, ok?
            self._granule[i] += self._deltaPhi[i] * dt / 100



'''
# average psi6 phases for each granule in frame 26, frame 27
# source: psi6_tracked_avg_plot_granule.m
# order of granules: big grain to the left, granules 1-10, big grain to the right
phis = [[-2.154, -2.154], 
        [-0.8657, -2.4625],
        [-1.6750, -2.4488],
        [-1.1998, -2.2778],
        [-1.0643, -2.2407],
        [-1.3656, -1.7055],
        [-0.7467, 1.0682],
        [0.1773, 1.4679],
        [-1.2103, 1.2222],
        [-1.0903, 1.0500],
        [-0.5810, 0.7767],
        [1.071, 1.071]]

# to convert psi6 phase to average misorientation angle, divide by 6
inits = [ phis[x][0]/6.0 for x in range(len(phis)) ]
finals = [ phis[x][1]/6.0 for x in range(len(phis)) ]

GBs = []
with open('gbs.csv') as csvFile:
    reader = csv.reader(csvFile, delimiter='\t')
    for row in reader:
        GBs.append( [ float(row[0]),int(row[1]),int(row[2])] )



numSteps = 200
time = [t for t in range(1+numSteps)]

myCrystal = Crystal(inits,finals,GBs)

energies = [myCrystal.totalEnergy()]
granEnergies = np.zeros((len(phis), 1+numSteps))
grangles = np.zeros((len(phis), 1+numSteps))

for gran in range(len(phis)):
    granEnergies[gran][0] = myCrystal.granuleEnergy(gran)
    grangles[gran][0] = myCrystal._granule[gran]

for i in range(numSteps):
    myCrystal.advanceInTime(100/numSteps)
    for gran in range(len(phis)):
        granEnergies[gran][i+1] = myCrystal.granuleEnergy(gran)
        grangles[gran][i+1] = myCrystal._granule[gran]
    energies.append(myCrystal.totalEnergy())
plt.plot(time, energies)
plt.figure()



granule_rgb_poster = [ (0.8863, 0.5961, 0.2706),
                (0.7725, 0.2392, 0.1686),
                (0.9255, 0.8118, 0.2784),
                (0.9216, 0.3804, 0.4980),
                (0.8471, 0.4196, 0.1647),
                (0.3647, 0.5216, 0.8275),
                (0.6196, 0.4275, 0.7333),
                (0.3020, 0.1725, 0.6157),
                (0.4314, 0.4471, 0.9176),
                (0.2549, 0.2314, 0.7961) ]

granule_rgb = plt.cm.nipy_spectral(np.linspace(0,1,10))  

labels = []
for gran in range(1,11):
    plt.plot(time, granEnergies[gran],color=granule_rgb[gran-1])
    labels.append(gran)
plt.legend(labels)
plt.figure()

for gran in range(1,11):
    plt.plot(time, grangles[gran],color=granule_rgb[gran-1])
plt.legend(labels)
plt.show()'''




plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 

a = 0 * np.pi/180
b = 40 * np.pi/180
phi0 = 20 * np.pi/180

inits = [a,phi0,phi0,b]
finals = [a,a,b,b]
GBs = [[0,1],[1,2],[2,3]]

numSteps = 200
time = [t/200*30 for t in range(1+numSteps)]

ds = [5, 10, 20, 50, 100]
ds = [50]
deltaE = []
for d in ds:
    GBs = [ [np.pi/2*d,0,1] , [d,1,2], [np.pi/2*d,2,3]]
    myCrystal = Crystal(inits,finals,GBs)
    energies = [myCrystal.totalEnergy()]
    for i in range(numSteps):
        myCrystal.advanceInTime(100/numSteps)
        energies.append(myCrystal.totalEnergy())
    deltaE.append(max(energies)-energies[0])
    maxi = 0
    for i in range(numSteps):
        if energies[i] == max(energies):
            maxi = i
    minE = min(energies)
    scaleMe = max(energies)-minE
    energies_shifted = [(x-minE)/scaleMe for x in energies]
    plt.plot(time, energies_shifted, color='green',linewidth=3)

    myXTix = [0,30,maxi/200*30]
    myYTix = [energies_shifted[0],max(energies_shifted),energies_shifted[-1]]

    myXTixRounded = [int(x*100)/100 for x in myXTix]
    myYTixRounded = [int(x*100)/100 for x in myYTix]

    

    plt.xticks(myXTixRounded)
    plt.yticks(myYTixRounded)

    
#plt.legend(ds)


#plt.figure()
#areas = [0.25*np.pi*d*d for d in ds]
#plt.plot(areas,deltaE)

plt.show()
