# 3 dimension fdtd acoustic simulation
from pylab import *
from scipy import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

X=40
Y=30

dx=0.1
dy=0.1
dt=0.0001

Ro=1.21
C=343
K=Ro*C*C

ang = arange(-pi,pi,2*pi/50)
sig = cos(ang)
sig += 1

P1 = zeros((X,Y),"float64")
P2 = zeros((X,Y),"float64")

Ux1 = zeros((X+1,Y),"float64")
Ux2 = zeros((X+1,Y),"float64")

Uy1 = zeros((X,Y+1),"float64")
Uy2 = zeros((X,Y+1),"float64")

#fig = figure()
#ax = axes3d.Axes3D(fig)
#
#x = arange(0, dx*X,dx)
#y = arange(0, dy*Y,dy)
#xx, yy = meshgrid(y,x)
#
#surf = ax.plot_surface(xx,yy,P2,rstride=1, cstride=1,cmap=cm.jet)
#oldsurf=surf

Z=Ro*C*20

for t in range(800):
    if t<len(sig):
        P1[X/2,Y/2] += sig[t]

    Ux2[1:X,:]=Ux1[1:X,:]-dt/Ro/dx*(P1[1:X,:]-P1[:X-1,:])
    Uy2[:,1:Y]=Uy1[:,1:Y]-dt/Ro/dy*(P1[:,1:Y]-P1[:,:Y-1])
    
    A=dt/(Ro*dx)
    #Ux2[0,:]=P1[0,:]/-Z
    #Ux2[-1,:]=P1[-1,:]/Z
    Ux2[0,:]=(1./(1.+Z*A))*( (1.-Z*A)*Ux1[0,:]-2*A*P1[0,:])
    Ux2[-1,:]=(1./(1.+Z*A))*( (1.-Z*A)*Ux1[-1,:]+2*A*P1[-1,:])

    
    #Uy2[:,0]=P1[:,0]/-Z
    #Uy2[:,-1]=P1[:,-1]/Z
    Uy2[:,0]=(1./(1.+Z*A))*( (1.-Z*A)*Uy1[:,0]-2*A*P1[:,0])
    Uy2[:,-1]=(1./(1.+Z*A))*( (1.-Z*A)*Uy1[:,-1]+2*A*P1[:,-1])

    P2[:X,:Y] = P1[:X,:Y]-K*dt/dx*(Ux2[1:X+1,:]-Ux2[:X,:]) \
                -K*dt/dy*(Uy2[:,1:Y+1]-Uy2[:,:Y])

#    surf = ax.plot_surface(xx,yy,P2,rstride=1, cstride=1,cmap=cm.jet)
#    ax.set_zlim3d(-1,1)
#    ax.collections.remove(oldsurf)
#    oldsurf=surf
#    draw()
#
    P1,P2=P2,P1
    Ux1,Ux2=Ux2,Ux1
    Uy1,Uy2=Uy2,Uy1

fig = figure()
ax = axes3d.Axes3D(fig)

x = arange(0, dx*X,dx)
y = arange(0, dy*Y,dy)
xx, yy = meshgrid(y,x)

surf = ax.plot_surface(xx,yy,P2,rstride=1, cstride=1,cmap=cm.jet)
oldsurf=surf
