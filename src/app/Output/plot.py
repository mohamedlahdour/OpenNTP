#!/home/tflaspoehler/python/bin/python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Arc
from matplotlib.patches import Circle, PathPatch
from matplotlib.transforms import Bbox
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
 
file_name = "pos20"
positions = np.loadtxt(file_name)
if (file_name[3:]=='2'):
    x = np.array( [positions[0]] )
    y = np.array( [positions[1]] )
    z = np.array( [positions[2]] )
    w = np.array( [positions[3]] )
else:
    x = positions[:,0]
    y = positions[:,1]
    z = positions[:,2]
    w = positions[:,3]
print "sum(w)",w.sum()
w = 50. * w / w.max()
     
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim3d([0, 1])
ax.set_ylim3d([0, 1])
ax.set_zlim3d([0, 1])
ax.set_autoscale_on(False)
n = 100
zr = []
for i in range(0, len(x)):
    add = True
    for j in range(0, len(zr)):
        if (zr[j][0] == z[i]):
            add = False
    if (add==True):
        zr.append( [z[i], ((x[i]**2.)+(y[i]**2.))**0.5] )
 
circle = []
box = Bbox(((-1,-1),(1,1)))
for i in range(0, len(zr)):
    r = zr[i][1]
    d = 2*r
    zed = zr[i][0]
    print "r,d,z",r,d,zed
    circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="Black", alpha=0.2))
    ax.add_patch(circle[-1])
    art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="x")
    circle[-1].set_clip_on(True)
    
    circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="Black", alpha=0.2))
    ax.add_patch(circle[-1])
    art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="y")
    circle[-1].set_clip_on(True)
    
    circle.append(Arc((0,0), width=d, height=d, angle=0.0, theta1=0.0, theta2=90.0, color="Black", alpha=0.2))
    ax.add_patch(circle[-1])
    art3d.pathpatch_2d_to_3d(circle[-1], z=zed, zdir="z")
    circle[-1].set_clip_on(True)
 
for i in range(0, len(x)):
    custom = ( x[i], y[i], z[i] )
    ax.scatter(x[i],y[i],z[i],color=custom,s=w[i]*5,zorder=2)
ax.view_init(elev=20., azim=25)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
 
plt.title(r'$S_{'+''.join([i for i in file_name[3:]])+'}$')
fig.savefig(file_name+".png", dpi=100, facecolor='white')
plt.show()
