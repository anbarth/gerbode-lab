import numpy as np

###    myGeo is where i keep the useful geometry functions i pillage from stackoverflow

# function which takes an origin point and returns {a function that gives a point's angle and distance from the origin}
# this is useful for sorting points in CW or CCW order
# stackoverflow.com/questions/41855695/sorting-list-of-two-dimensional-coordinates-by-clockwise-angle-using-python
def make_clockwiseangle_and_distance(origin):
    def clockwiseangle_and_distance(point):
        refvec = [0,1]
        # Vector between point and the origin: v = p - o
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
    return clockwiseangle_and_distance

# shoelace formula for calculating the area of a polygon
# stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
def polyArea(corners):
    n = len(corners) # of corners
    area = 0.0
    for i in range(n):
        j = (i + 1) % n
        area += corners[i][0] * corners[j][1]
        area -= corners[j][0] * corners[i][1]
    area = abs(area) / 2.0
    return area

# find the crossing points of two circles
# one circle of radius r0 at (x0,y0) and another of radius r1 at (x1,y1)
# stackoverflow.com/questions/55816902/finding-the-intersection-of-two-circles
# if no crossing points, returns None
# if 2 crossing points, returns both in the form (x1,y1,x2,y2)
# if the circles just kiss, returns the singular crossing point twice (x1,y1,x1,y1)
def circIntersections(x0, y0, r0, x1, y1, r1):
    # circle 1: (x0, y0), radius r0
    # circle 2: (x1, y1), radius r1

    d=np.sqrt((x1-x0)**2 + (y1-y0)**2)
    
    # non intersecting
    if d > r0 + r1 :
        return None
    # One circle within other
    if d < abs(r0-r1):
        return None
    # coincident circles
    if d == 0 and r0 == r1:
        return None
    else:
        a=(r0**2-r1**2+d**2)/(2*d)
        h=np.sqrt(r0**2-a**2)
        x2=x0+a*(x1-x0)/d   
        y2=y0+a*(y1-y0)/d   
        x3=x2+h*(y1-y0)/d     
        y3=y2-h*(x1-x0)/d 

        x4=x2-h*(y1-y0)/d
        y4=y2+h*(x1-x0)/d
        
        return (x3, y3, x4, y4)

# finds the centroid of a set of points
def centroid(x_coords,y_coords):
    _len = len(x_coords)
    centroid_x = sum(x_coords)/_len
    centroid_y = sum(y_coords)/_len
    return [centroid_x, centroid_y]