import pylab as pl
import mpl_toolkits.mplot3d.axes3d as p3

#data is an ndarray with the necessary data and colors is an ndarray with
#'b', 'g' and 'r' to paint each point according to its class
fig=pl.figure()
ax = p3.Axes3D(fig)
data=pl.loadtxt('p.dat')
i=0
x=[]
y=[]
z=[]
for pos in data[:]: 
  
  if ((i/3)*3==i):
    x=x+[pos[0]]
    y=y+[pos[1]]
    z=z+[pos[2]]
  i=i+1
ax.scatter(x,y,z,s=5,color='red')

data=pl.loadtxt('p-all.dat')
i=0
x=[]
y=[]
z=[]
for pos in data[:]: 
  
  if ((i/3)*3==i):
    x=x+[pos[0]]
    y=y+[pos[1]]
    z=z+[pos[2]]
  i=i+1

ax.scatter(x,y,z,s=1,c='yellow')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
fig.add_axes(ax)
pl.show()
