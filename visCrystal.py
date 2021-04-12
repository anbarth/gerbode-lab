
import csv
import matplotlib.pyplot as plt

fname = 'crystals_oldformat/annasNewCrystal.csv'

particleCenters = []
with open(fname) as csvFile:
    reader = csv.reader(csvFile, delimiter=',')
    first = True
    for row in reader:
        # first line of csv gives xmin, ymin, window dimensions, and bead radius
        if first:
            beadRad = float(row[4])
            xmin = float(row[0])
            ymin = float(row[1])
            windowSize = [float(row[2]),float(row[3])]
            first = False
            continue
        # after the first line, it's particle centers all the way down
        p = ( (float(row[0])) , 
                (float(row[1])) )
        particleCenters.append(p)


fig, ax = plt.subplots()
ax.set_aspect(1)

for p in particleCenters:
    circ = plt.Circle(p, beadRad, facecolor='gray',edgecolor=None, alpha=0.3)
    ax.add_artist(circ)
    circ = plt.Circle(p, beadRad/20, facecolor='k', edgecolor=None)
    ax.add_artist(circ)


plt.xlim(xmin,xmin+windowSize[0])
plt.ylim(ymin,ymin+windowSize[1])

plt.show()