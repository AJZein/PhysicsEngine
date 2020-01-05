import numpy as np

def slope(point1, point2):
    return (point2[1]-point1[1])/(point2[0]-point1[0])

def dist(p1, p2):
    """Returns distance between two vectors of length 2"""
    value = (p1[0] - p2[0])**2
    value += (p1[1] - p2[1])**2
    return (value)**(0.5)

def norm(vector):
    return (vector[0]**2 + vector[1]**2)**(0.5)

def drawPoly(coords, line):
    """Draws a polygon, given a matplotlib line object"""
    x = coords[:,0]
    y = coords[:,1]
    x = np.append(x, x[0])
    y = np.append(y, y[0])
    line.set_data(x, y)
    return line