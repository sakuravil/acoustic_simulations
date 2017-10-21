from time import *
from pylab import *
from scipy import *
from numpy import *
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.pyplot as plt
import csv

X=101
Y=101
NT=301

dx=5.e-2
dy=5.e-2
dt=5.e-5

Ro=1.21
C=343
K=Ro*C*C

P1 = zeros((X,Y),"float")
Ux = zeros((X+1,Y),"float")
Uy = zeros((X,Y+1),"float")

sigma = 0.2
xc= (X-1) // 2
yc= (Y-1) // 2
for i in range(1, X-1):
    x = dx * i
    for j in range(1, Y-1):
        y = dy * j
        TX = x - xc * dx
        TY = y - yc * dy
        if j <= 60 and j >= 40:
            P1[i][j] = 1
        else:
            P1[i][j] = 0

Z=Ro*C

start=clock()

list=[]

for t in range(NT):
    if(t%50)==0:
        for var in range(X):
            list.append(P1[xc][var])
        f = open('fdtd.csv', 'a')
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(list)
        f.close()
        list = []

    Ux[1:X,:]=Ux[1:X,:]-dt/Ro/dx*(P1[1:X,:]-P1[:X-1,:])
    Uy[:,1:Y]=Uy[:,1:Y]-dt/Ro/dy*(P1[:,1:Y]-P1[:,:Y-1])
    
    Ux[0,:]=P1[0,:]/-Z
    Ux[-1,:]=P1[-1,:]/Z
    
    Uy[:,0]=P1[:,0]/-Z
    Uy[:,-1]=P1[:,-1]/Z

    P1[:X,:Y] = P1[:X,:Y]-K*dt/dx*(Ux[1:X+1,:]-Ux[:X,:]) \
                -K*dt/dy*(Uy[:,1:Y+1]-Uy[:,:Y])

end=clock()

print ("Processing Time : " + str((end - start)) + " [sec]")
