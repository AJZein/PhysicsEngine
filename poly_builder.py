import numpy as np
from tools import dist

class Polygon:
    """Class for storing all polygon data"""
    
    def __init__(self, points):
        """Points is an nx2 numpy array containing all the points which, when joined sequentially, form the polygon"""
        self.sides = points.shape[0]
        
        self.vertices = np.zeros_like(points) # The vertices will be relative to the centre of geometry
        self.centre = np.sum(points, axis=0)/self.sides
        for i in range(self.sides):
            self.vertices[i] = points[i] - self.centre # The vertices are relative to the centre of geometry
        
        self.mass = 1
        self.force = np.zeros((2))
        self.acceleration = np.zeros((2))
        self.velocity = np.zeros((2))
        self.angVelocity = 0
        self.radius = max([dist([0,0], vertex) for vertex in self.vertices])
        
        
    def update(self, dt):
        self.acceleration = self.force/self.mass
        self.velocity += dt*self.acceleration
        self.centre += dt*self.velocity
        # Angular velocity limiter
        #if abs(self.angVelocity) > 1:
        #    self.angVelocity /= abs(self.angVelocity)
        self.rotate(dt*self.angVelocity)
        
    def getPos(self):
        """Returns the absoulute positions of all points in the polygon"""
        positions = np.zeros_like(self.vertices)
        for i in range(self.sides):
            positions[i] = self.vertices[i] + self.centre
        return positions
    
    def getPosCirc(self):
        """Returns the absoulute positions of all points in the polygon including the first point appended to the end"""
        positions = self.getPos()
        positions = np.vstack((positions, positions[0]))
        return positions
    
    def rotate(self, angle, degree=False, centre=False):
        """Perform counter-clockwise rotation of a set of coords about that
    coords centre of geometry"""
    
        if degree:
            angle = 2*np.pi*angle/360
            
        Rmatrix = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

        if not centre:
            centre = np.sum(self.vertices, axis=0)/self.sides # Centre of geometry for vertices
            
        rotCentre = centre - np.matmul(Rmatrix, centre)
        for i in range(self.sides):
            self.vertices[i] = np.matmul(Rmatrix, self.vertices[i]) + rotCentre


def createPoints(lower,upper,decimals, number=1):
    if number == 1:
        return np.around(np.random.rand(2)*(upper-lower) + lower, decimals=decimals)
    else:
        return np.around(np.random.rand(number, 2)*(upper-lower) + lower, decimals=decimals)


def checkInter(points, decimals=3):
    """Checks if the line formed by the first two points intersects the line formed by the other two points. `points` is a 4x2 numpy array [[x1, y1], [x2, y2]...]"""
    
    if points[1, 0]-points[0, 0]:
        m1 = (points[1, 1] - points[0, 1])/(points[1, 0]-points[0, 0])
        c1 = points[0, 1] - m1*points[0, 0]
    else:
        m1 = 0
        c1 = points[1, 1] - points[0, 1]
        
    if points[3, 0]-points[2, 0]:
        m2 = (points[3, 1] - points[2, 1])/(points[3, 0]-points[2, 0])
        c2 = points[2, 1] - m2*points[2, 0]
    else:
        m2 = 0
        c2 = points[1, 1] - points[0, 1]
    
    if m1 == m2:
        return False
    
    x = round((c2 - c1)/(m1 - m2), decimals)
    
    if min(points[0,0], points[1,0]) < x < max(points[0,0], points[1,0]):
        if min(points[2,0], points[3,0]) < x < max(points[2,0], points[3,0]):
            return True
    return False

def checkInterLine(line,array, simple=True):
    """Check if a line intersects with any of the lines formed by the connecting the points in array sequentially"""
    
    numLines = array.shape[0]
    array = np.vstack((array, array[0]))
    points = np.zeros((4, 2))
    points[0:2] = line

    if not simple:
        intersectingPoints = []

    for i in range(numLines):
        if not np.array_equal(points[0], array[i]):
            points[2] = array[i]
            points[3] = array[i+1]
            if checkInter(points):
                if simple:
                    return True
                else:
                    intersectingPoints.append(points)

    if simple:
        return False
    else:
        return intersectingPoints

def checkInterSet(array1, array2):
    """Checks if the polygon formed by connecting all the points in array1 sequentially intersects a polygon generated from array2 in the same way"""
    
    for i in range(array1.shape[0]-1):
        if checkInterLine(array1[i:i+2], array2):
            return True
        
    if checkInterLine(np.vstack((array1[-1], array1[0])), array2):
        return True
        
    return False

def findInterSet(array1, array2):
    """Returns a list of 4x2 numpy arrays which represent all the intersecting lines between array1 and array2. Each row contains a point. 
    The first two rows in each array form a line from array1 and the final two rows form a line from array2"""
    
    allIntersections = []
    for i in range(array1.shape[0]-1):
        intersects = checkInterLine(array1[i:i+2], array2, simple=False)
        if intersects:
            allIntersections.extend(intersects)

    return allIntersections

def sortbyslope(array, origin=[0,0]):
    """Sorts a set of points according to their slopes from a given origin"""

    quadrants = {'quad1':[], 'quad2':[], 'quad3':[], 'quad4':[]}
    axis = {'y+':[], 'y-':[]}
    
    ox, oy = origin[0], origin[1]
    for i in range(array.shape[0]):
        point = array[i]
        if point[0] != ox:
            slope = (point[1] - oy)/(point[0] - ox)
        else:
            if point[1] > 0:
                axis['y+'].append(point) # Points on y+ are sorted here
            else:
                axis['y-'].append(point) # Points on y- are sorted here
            
        if point[0] > ox and point[1] >= oy: # Points on x+ are sorted here
            quadrants['quad1'].append((slope, point))
        elif point[0] < ox and point[1] >= oy: # Points on x- are sorted here
            quadrants['quad2'].append((slope, point))
        elif point[0] < ox and point[1] < oy:
            quadrants['quad3'].append((slope, point))
        elif point[0] > ox and point[1] < oy:
            quadrants['quad4'].append((slope, point))
        
    sortedArray = np.zeros_like(array)
    n = 0
    for i in range(1, 5):
        quadrants['quad'+str(i)].sort(key=lambda x:x[0])
        for j in range(len(quadrants['quad'+str(i)])):
            sortedArray[j+n] = quadrants['quad'+str(i)][j][1]
        n += len(quadrants['quad'+str(i)])
        
        if i == 1 and axis['y+']:
            for j in range(len(axis['y+'])):
                sortedArray[j+n] = axis['y+'][j]
            n += len(axis['y+'])
            
        if i == 3 and axis['y-']:
            for j in range(len(axis['y-'])):
                sortedArray[j+n] = axis['y-'][j]
            n += len(axis['y-'])
        
    return sortedArray
        
def rand_polygon(sides, lower=0, upper=5, mutate=False):
    """Creates a random polygon with a given number of sides which lies entirely between `lower` and `upper` in both the x and y axis"""
    origin = [(upper+lower)/2, (upper+lower)/2]
    points = createPoints(lower, upper, 3, number=sides)
    points = sortbyslope(points, origin = origin)
    n = 0
    while checkInterSet(points, points): # Check if the polygon intersects itself
        points = createPoints(lower, upper, 3, number=sides)
        points = sortbyslope(points, origin = origin)
        n += 1
        if n >= 10:
            raise Exception("Error: Too many invalid polygon generated")
    
    if mutate:
        points = mutation(points, lower, upper)
        
    return points
        
def mutation(points, lower, upper, rounds=20, sticky=0.2):
    """Mutates a simple polygon into another simple polygon"""
    
    numPoints = points.shape[0]
    mutations = createPoints(lower, upper/10, 3, number=points.shape[0])
    for run in range(rounds):
        
        for i in range(1, numPoints-1):
            points[i] += mutations[i]
            if points[i, 0] < lower or points[i, 0] > upper or points[i, 1] < lower or points[i, 1] > upper:
                points[i] -= mutations[i]
            elif checkInterLine(points[i:i+2], points) or checkInterLine(points[i-1:i+1], points):
                points[i] -= mutations[i]
            else:
                mutations[i+1] += sticky*mutations[i]
                
        for i in range(2):
            points[i - 1] += mutations[i - 1]
            if points[i-1, 0] < lower or points[i-1, 0] > upper or points[i-1, 1] < lower or points[i-1, 1] > upper:
                points[i-1] -= mutations[i-1]
            elif checkInterLine(np.vstack((points[i - 1], points[i])), points) or checkInterLine(np.vstack((points[i - 2], points[i-1])), points):
                points[i - 1] -= mutations[i - 1]

    return points
