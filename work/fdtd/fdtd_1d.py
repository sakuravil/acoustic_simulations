# 1 dimension fdtd simulation
from pylab import *
from scipyq import *
from matplotlib import cm

X=40
Y=30

dx=1
dy=1
dt=0.0001

Ro=1.21
C=343
K=Ro*C*C

Q = 0.5+0.5*cos(arange(-pi,pi,2*pi/200))

P1 = zeros((X,Y),"float64")
P2 = zeros((X,Y),"float64")

Ux1 = zeros((X+1,Y),"float64")
Ux2 = zeros((X+1,Y),"float64")

Uy1 = zeros((X,Y+1),"float64")
Uy2 = zeros((X,Y+1),"float64")

mic=[]

for n in range(5000):
    if n<len(Q):
        P1[20,15] += Q[n]

    mic.append(P1[0,0])

    for x in range(X-1):
        for y in range(Y):
           Ux2[x+1,y]=Ux1[x+1,y]-dt/Ro/dx*(P1[x+1,y]-P1[x,y])
           
    for x in range(X):
        for y in range(Y-1):
            Uy2[x,y+1]=Uy1[x,y+1]-dt/Ro/dy*(P1[x,y+1]-P1[x,y])

    for x in range(X):
        for y in range(Y):
            P2[x,y] = P1[x,y]-K*dt/dx*(Ux2[x+1,y]-Ux2[x,y]) \
                -K*dt/dy*(Uy2[x,y+1]-Uy2[x,y])

    P1,P2=P2,P1
    Ux1,Ux2=Ux2,Ux1
    Uy1,Uy2=Uy2,Uy1

figure()
plot(mic)
xlabel("Sample")
ylabel("Relative sound pressure")
show()

figure()
contourf(P1.T,aspect="equal", cmap=cm.jet)
xlim(0,X-1)
ylim(0,Y-1)
xlabel("X [sample]")
ylabel("Y [sample]")
show()
