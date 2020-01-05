import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from poly_builder import rand_polygon, Polygon
from collisions import checkBorders, checkCollisions
from tools import drawPoly
import random
  
def draw(polygons, energyLoss=0, save=False):
    """Creates an animation for a list of polygon objects"""
    dt = 0.02
    fig = plt.figure(0, figsize=(6,6))
    ax = plt.axes(xlim=(borderxm, borderxp), ylim=(borderym, borderyp))
    lines = []
    for i in range(len(polygons)):
        lineObj, = ax.plot([], [], lw=2)
        lines.append(lineObj)
        
    def init():
        return lines

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=120*4, interval=20, blit=True, 
                                   fargs=(polygons, dt, lines, energyLoss))
    

    if save:
        plt.rcParams['animation.ffmpeg_path'] = "C:\\Users\\alize\Programs\\ffmpeg-20191201-637742b-win64-static\\bin\\ffmpeg.exe"
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800, extra_args=['-vcodec', 'libx264'])
        anim.save('multi_animation.mp4', writer=writer, )

    plt.show()
        
def animate(frame, polygons, dt, lines, energyLoss):
    """Function used to draw each frame of animation"""
    for i in range(4):
        polygons = checkCollisions(polygons)
        polygons = checkBorders(polygons, [borderxm, borderxp, borderym, borderyp], energyLoss)
        for polygon in polygons:
            polygon.update(dt)

        
    for i in range(len(polygons)):
        lines[i] = drawPoly(polygons[i].getPos(), lines[i])
        
    return lines

nsides = 7
nPolys = 7
lower, upper = 0, 5
energyLoss = 0.0 # Number between 0 and 1 for amount of velocity lost in each collision
borderxm, borderxp, borderym, borderyp = 0, 20, 0, 20



polygons = []
gridPoints = int((nPolys)**(0.5))+3
gridx, gridy = (borderxp - borderxm)/(gridPoints), (borderyp - borderym)/(gridPoints)
n = 0
for i in range(1, gridPoints-1):
    for j in range(1, gridPoints-1):
        sides = nsides + random.randint(-3, 3)
        polygon = Polygon(rand_polygon(sides, lower=lower, upper=upper))
        polygon.centre = [gridx*j, gridy*i]
        polygon.velocity = [random.randint(-3, 3), random.randint(-3, 3)]
        polygons.append(polygon)
        n += 1
        if n >= nPolys:
            break

borders = Polygon(np.array([[borderxm,borderym], [borderxm, borderyp], 
                            [borderxp,borderyp], [borderxp,borderym]]))

draw(polygons, energyLoss=energyLoss)