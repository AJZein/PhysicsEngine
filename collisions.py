import numpy as np
from poly_builder import findInterSet, checkInterLine
from tools import slope, dist, norm
    
def assignQuadrant(origin, point):
    """Assigns a quadrant to `point` using `origin` as origin"""
    point = point - origin
    if point[0] > 0 and point[1] >= 0:
        return 0
    elif point[0] <= 0 and point[1] > 0:
        return 1
    elif point[0] < 0 and point[1] <= 0:
        return 2
    elif point[0] >= 0 and point[1] < 0:
        return 3
    else:
        raise Exception(f"Something strange happened, origin was {origin}, point was {point}")

def isInside(point, vertices):
    """Checks if `point` is inside the polygon formed by `vertices`"""
    winding = 0

    for i in range(vertices.shape[0] + 1):
        if i != vertices.shape[0]:
            vertex = vertices[i]
        else:
            vertex = vertices[0] # Connect last point to initial point
            
        if not np.any(vertex - point): # Check whether the point is one of our vertices
            return False
        
        newQuadrant = assignQuadrant(point, vertex)
        # Increasing winding number
        if i > 0:
            if newQuadrant == (oldQuadrant + 1)%4: # Counter-clockwise windings increase winding number positvely
                winding += 0.25
            elif newQuadrant == (oldQuadrant - 1)%4:
                winding -= 0.25
            elif newQuadrant != oldQuadrant: 
                # If the polygon has connections that span 2 quadrants, extra calculations are needed to 
                # determine whether the winding is clockwise or counterclockwise
                oldVertex = vertices[i-1]
                if (oldVertex[0] == point[0] and oldVertex[0] == vertex[0]): # Checking if the vertices and point form a vertical line
                    return False
                elif (oldVertex[0] == point[0]): # Checking if the old vertex froms a vertical line with the point
                    winding += 0.5
                    continue


                # Slopes are computed now since there is no chance of an infinite slope and used to determine winding direction
                m1, m2 = slope(oldVertex, point), slope(oldVertex, vertex)
                if m1 > m2:
                    winding += 0.5
                elif m1 < m2:
                    winding -= 0.5
                else: # This implies the point lies on an edge
                    return False
                    

        oldQuadrant = newQuadrant

    if abs(winding) >= 1:
        return True
    
    return False

def isIntersecting(polygons):
    """Checks whether two polygons are intersecting and returns the indices of the points which are intersecting in each"""
    #Setup
    positionSets = []
    for polygon in polygons:
        positionSets.append(polygon.getPos())

    #Finding which points are inside
    intersectionSets = [[], []]
    for i in range(len(positionSets[0])):
        point = positionSets[0][i]
        if isInside(point, positionSets[1]):
            intersectionSets[0].append(i)

    for i in range(len(positionSets[1])):
        point = positionSets[1][i]
        if isInside(point, positionSets[0]):
            intersectionSets[1].append(i)

    return intersectionSets

def flipPoints(polygons, intersectionSets):
    """Computes the direction of the flip vectors based on the points which lie inside each polygon"""
    # Computing velocity flip direction
    flipVectors = []
    for polygon, indices in zip(polygons, intersectionSets):
        if indices:
            line1 = polygon.getPosCirc()[indices[0]] - polygon.getPosCirc()[indices[0] - 1] # First line which intersect with the other polygon
            line2 = polygon.getPosCirc()[indices[-1]] - polygon.getPosCirc()[(indices[-1] + 1)] # Last line which intersects with the other polygon
                
            # Translational motion will be flipped along the average of the first and last intersecting line with the other polygon
            flipVector = (line1 + line2)/2 # This vector currently point towards the other polygon
            flipVectors.append(flipVector)
        else:
            flipVectors.append(0)

    if not intersectionSets[0]:
        flipVectors[0] = -1*flipVectors[1]
    elif not intersectionSets[1]:
        flipVectors[1] = -1*flipVectors[0]

    return flipVectors

def flipPolys(polygons):
    """Returns two lists, each containing vectors along which the respectives polygons must have their velocities flipped"""
    intersections = findInterSet(polygons[0].getPos(), polygons[1].getPos())
    flipVectors = [[], []]

    if not intersections:
        return flipVectors

    for points in intersections:
        for i in range(2):
            
            mOrg = (points[1 + i*2, 1] - points[0 + i*2, 1])/(points[1 + i*2, 0]-points[0 + i*2, 0])
            mPerp = -1/mOrg
            cPerp = points[0 + i*2, 1] - mPerp*points[0 + i*2, 0]
            # Determining direction of new vector
            yPerp1 = mPerp*points[0 + i*2,0] + cPerp
            yPerp2 = mPerp*points[1 + i*2,0] + cPerp

            flipPoints = np.array([[points[0 + i*2,0], yPerp1], [points[1 + i*2,0], yPerp2]])

            if dist(flipPoints[0], polygons[i].centre) > dist(flipPoints[1], polygons[i].centre):
                flipVector = flipPoints[0] - flipPoints[1]
            else:
                flipVector = flipPoints[0] - flipPoints[1]
            flipVectors[i].append(flipVector)

    return flipVectors


def projection(basis1, basis2, vector):
    """Computes `vector` as a linear combination of two other basis vectors
    Vector = A*basis1 + B*basis2"""
    
    r1, r2 = vector[1]/basis1[1], basis2[1]/basis1[1]
    B = (vector[0] - r1)/(basis2[0] - basis1[0]*r2)
    A = r1 - B*r2
    
    return A, B

def projection_single(basis1, vector):
    """Computes `vector` as a linear combination of two other basis vectors. basis1 is provided and basis2 is calculated
    by rotating basis1 90 degrees. vector = A*basis1 + B*basis2"""

    R90 = np.array([[0, -1], [1, 0]])
    basis2 = np.matmul(R90, basis1)
    A, B = projection(basis1, basis2, vector)
    return A, B, basis2

def rebound(masses, velocities):
    """Computes recoil velocities of two objects colliding, conserving momentum and kinetic energy. Collision totally elastic"""
    m1, m2, v1, v2 = masses[0], masses[1], velocities[0], velocities[1]
    v11 = m1*v1 - m2*v1 + 2*m2*v2
    v11 /= (m1+m2)
    v22 = 2*m1*v1 - m1*v2 + m2*v2
    v22 /= (m1+m2)
    
    #print(m1*v1 + m2*v2 - m1*v11 - m2*v22) # Conservation of momentum
    #print(m1*(v1**2) + m2*(v2**2) - m1*(v11**2) - m2*(v22**2)) # Conservation of energy
    return [v11, v22]
    

def collision(polygons):
    """Adjusts velocity and rotation of two polygon objects if they collide"""

    intersectionSets = isIntersecting(polygons)
    if not intersectionSets[0] and not intersectionSets[1]:
        return

    # Computing velocity flip direction
    flipVectors = flipPoints(polygons, intersectionSets)
        
    # Alternate method of computing rebound direction, doesnt work well
    if False:
        flipVectors = flipPolys(polygons)
        if flipVectors[0]:
            return

        flipVectors = [sum(flipVectors[0])/len(flipVectors[0]), sum(flipVectors[1])/len(flipVectors[1])] # Repalce with some sort of averaging rather than taking just the first entry

    # Flipping Velocity
    newVelocities = []
    for basis1, polygon in zip(flipVectors, polygons):
        flipA, flipB, basis2 = projection_single(basis1, polygon.velocity + polygon.angVelocity)
        newVelocities.append(-abs(flipA)*basis1 + flipB*basis2)

    # Normalising velocity to conserve momentum and kinetic energy
    oldVelMags = [norm(polygons[i].velocity) for i in range(2)] # Magnitude of velocities before collision
    newVelMags = [norm(newVelocities[i]) for i in range(2)] # When flipped the magnitude of the vectors change to these values
    reboundVelMags = rebound([polygons[0].mass, polygons[1].mass], oldVelMags) # New vectors should have these magnitudes to conserve momentum
    
    if newVelMags[0] == 0:
        newVelocities[0] = 0
    else:
        newVelocities[0] = newVelocities[0]*reboundVelMags[0]/newVelMags[0]
    if newVelMags[1] == 0:
        newVelocities[1] = 0
    else:
        newVelocities[1] = newVelocities[1]*reboundVelMags[1]/newVelMags[1]

    #normConstants = [reboundVelMags[i]/newVelMags[i] for i in range(2)] # Values to use to normalise the flipped vectors to the correct magnitudes
    #newVelocities = [newVelocities[i]*normConstants[i] for i in range(2)]

    # Computing rotation
    torques = []
    for i in range(2):
        torques.append(collisionTorque(polygons[i], newVelocities[i], intersectionSets[i]))

    if not intersectionSets[0]:
        torques[0] = -1*torques[1]
    elif not intersectionSets[1]:
        torques[1] = -1*torques[0]
    
    # Setting new velocity and rotation
    for i in range(2):
        polygons[i].velocity = newVelocities[i]
        polygons[i].angVelocity = torques[i]
     
def collisionTorque(polygon, newVelocity, collisionIndices, fixedPoint=0):
    """"Computes torque for a polygon whose velocity is changed to `newVelocity` through a collision.
    Setting fixed point allows the collisionIndices to be set to 0 and for torque to be calculated about these points"""

    forces = np.ones_like(polygon.vertices)*polygon.velocity
    for i in collisionIndices:
        if isinstance(fixedPoint, np.ndarray):
            forces[i] *= 0
        else:
            forces[i] = newVelocity

    if isinstance(fixedPoint, np.ndarray):
        torque = np.sum(np.cross(polygon.vertices + (polygon.centre - fixedPoint), forces))
    else:
        torque = np.sum(np.cross(polygon.vertices, forces))
    return torque*0.2

def checkCollisions(polygons):
    """Checks which polygons might be colliding by computing centre to centre distance"""
    for i in range(len(polygons)-1):
        for j in range(i+1, len(polygons)):
            if dist(polygons[i].centre, polygons[j].centre) < polygons[i].radius + polygons[j].radius:
                collision([polygons[i], polygons[j]])         
    return polygons

def checkBorders(polygons, borders, energyLoss):
    """Collision algorithm for borders"""
    
    for polygon in polygons:
        positions = polygon.getPos()
        collisionPoints = []
        xdir, ydir = 0, 0
        for i in range(polygon.sides):
            if positions[i, 0] <= borders[0]:
                collisionPoints.append(i)
                xdir = 1
            elif positions[i, 0] >= borders[1]:
                collisionPoints.append(i)
                xdir = -1
                
            if positions[i, 1] <= borders[2]:
                collisionPoints.append(i)
                ydir = 1
            elif positions[i, 1] >= borders[3]:
                collisionPoints.append(i)
                ydir = -1

        if collisionPoints:
            newVelocity = polygon.velocity
            if xdir:
                newVelocity[0] = xdir*abs(newVelocity[0])

            if ydir:
                newVelocity[1] = ydir*abs(newVelocity[1])

            polygon.velocity = newVelocity*(1-energyLoss)
            torque = collisionTorque(polygon, newVelocity, collisionPoints, positions[collisionPoints[0]])
            polygon.angVelocity = torque/10
         


    return polygons