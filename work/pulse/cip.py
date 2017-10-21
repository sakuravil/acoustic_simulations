from time  import *
from pylab import *
from numpy import *
from scipy import *
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import csv

X = 101
Y = 101
NT = 301

f_p   = zeros((X,Y), "float64")
f_m   = zeros((X,Y), "float64")
g_p   = zeros((X,Y), "float64")
g_m   = zeros((X,Y), "float64")
fn_p  = zeros((X,Y), "float64")
fn_m  = zeros((X,Y), "float64")
gn_p  = zeros((X,Y), "float64")
gn_m  = zeros((X,Y), "float64")

P     = zeros((X,Y), "float64")
dx_P  = zeros((X,Y), "float64")
dy_P  = zeros((X,Y), "float64")

ZUx    = zeros((X,Y), "float64")
dx_ZUx = zeros((X,Y), "float64")

ZUy    = zeros((X,Y), "float64")
dy_ZUy = zeros((X,Y), "float64")

coeff10 = zeros(8,  "float64")
coeff31 = zeros(32, "float64")

xc= (X-1) // 2
yc= (Y-1) // 2

dx = 5.e-2
dy = 5.e-2
dt = 5.e-5

Ro = 1.21
bm = 1.4235529e5
c0 = sqrt(bm / Ro)
Z0 = sqrt(bm * Ro)
sigma = 0.2

Ua = c0
xi =-Ua * dt
C  = c0 * dt / dx
C2 = C  * C
C3 = C2 * C

coeff10[0]  = C * 0.5
coeff10[1]  = (-C + 1.) * 0.5

coeff31[0]  = (-2. * C3 + 3. * C2) * 0.5
coeff31[1]  = (2. * C3 - 3. * C2 + 1.) * 0.5
coeff31[2]  = xi * (C2 - C) * 0.5
coeff31[3]  = xi * (C2 - 2. * C + 1.) * 0.5
coeff31[4]  = 6. * (-C3 + C2) / xi * 0.5
coeff31[5]  = 6. * (C3 - C2) / xi * 0.5
coeff31[6]  = (3. * C2 - 2. * C) * 0.5
coeff31[7]  = (3. * C2 - 4. * C + 1.) * 0.5

Ua =-c0
xi =-Ua * dt
C  = c0 * dt / dx
C2 = C  * C
C3 = C2 * C

coeff10[2]  = C * 0.5
coeff10[3]  = (-C + 1.) * 0.5

coeff31[8]  = (-2. * C3 + 3. * C2) * 0.5
coeff31[9]  = (2. * C3 - 3. * C2 + 1.) * 0.5
coeff31[10] = xi * (C2 - C) * 0.5
coeff31[11] = xi * (C2 - 2. * C + 1.) * 0.5
coeff31[12] = 6. * (-C3 + C2) / xi * 0.5
coeff31[13] = 6. * (C3 - C2) / xi * 0.5
coeff31[14] = (3. * C2 - 2. * C) * 0.5
coeff31[15] = (3. * C2 - 4. * C + 1.) * 0.5

Ua = c0
yi =-Ua * dt
C  = c0 * dt / dy
C2 = C  * C
C3 = C2 * C

coeff10[4]  = C * 0.5
coeff10[5]  = (-C + 1.) * 0.5

coeff31[16] = (-2. * C3 + 3. * C2) * 0.5
coeff31[17] = (2. * C3 - 3. * C2 + 1.) * 0.5
coeff31[18] = yi * (C2 - C) * 0.5
coeff31[19] = yi * (C2 - 2. * C + 1.) * 0.5
coeff31[20] = 6. * (-C3 + C2) / yi * 0.5
coeff31[21] = 6. * (C3 - C2) / yi * 0.5
coeff31[22] = (3. * C2 - 2. * C) * 0.5
coeff31[23] = (3. * C2 - 4. * C + 1.) * 0.5

Ua =-c0
yi =-Ua * dt
C  = c0 * dt / dy
C2 = C  * C
C3 = C2 * C

coeff10[6]  = C * 0.5
coeff10[7]  = (-C + 1.) * 0.5

coeff31[24] = (-2. * C3 + 3. * C2) * 0.5
coeff31[25] = (2. * C3 - 3. * C2 + 1.) * 0.5
coeff31[26] = yi * (C2 - C) * 0.5
coeff31[27] = yi * (C2 - 2. * C + 1.) * 0.5
coeff31[28] = 6. * (-C3 + C2) / yi * 0.5
coeff31[29] = 6. * (C3 - C2) / yi * 0.5
coeff31[30] = (3. * C2 - 2. * C) * 0.5
coeff31[31] = (3. * C2 - 4. * C + 1.) * 0.5

for i in range(1, X-1):
    x = dx * i
    for j in range(1, Y-1):
        # y = dy * j
        # TX = x - xc * dx
        # TY = y - yc * dy
        if j <= 60 and j >= 40:
            P[i][j] = 1
        else:
            P[i][j] = 0
        dx_P[i][j]  = 0
        dy_P[i][j]  = 0
                    
def LINEAR(coeff0, coeff1, f0, f1):
    return    coeff0 * f0 \
            + coeff1 * f1

def CIP(coeff0, coeff1, coeff2, coeff3, f0, f1, g0, g1):
    return    coeff0 * f0 \
            + coeff1 * f1 \
            + coeff2 * g0 \
            + coeff3 * g1

start = clock()


list = []

for t in range(NT):
    if (t % 50) == 0:
        for var in range(X):
            list.append(P[xc][var])
        f = open('cip.csv', 'a')
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(list)
        f.close()
        list = []

    f_p[1:X-1,1:Y-1]   = dy_P[1:X-1,1:Y-1] + (ZUx[1:X-1,2:Y] \
                          - ZUx[1:X-1,0:Y-2]) / (2. * dy)
    f_m[1:X-1,1:Y-1]   = dy_P[1:X-1,1:Y-1] - (ZUx[1:X-1,2:Y] \
                          - ZUx[1:X-1,0:Y-2]) / (2. * dy)

    fn_p[1:X-1,1:Y-1]  = LINEAR(coeff10[0], coeff10[1], \
                                 f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]  = LINEAR(coeff10[2], coeff10[3], \
                                 f_m[2:X,1:Y-1],   f_m[1:X-1,1:Y-1])

    dy_P[1:X-1,1:Y-1]  = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) 

    f_p[1:X-1,1:Y-1]   =    P[1:X-1,1:Y-1] +ZUx[1:X-1,1:Y-1]
    f_m[1:X-1,1:Y-1]   =    P[1:X-1,1:Y-1] - ZUx[1:X-1,1:Y-1]
    g_p[1:X-1,1:Y-1]   = dx_P[1:X-1,1:Y-1] + dx_ZUx[1:X-1,1:Y-1]
    g_m[1:X-1,1:Y-1]   = dx_P[1:X-1,1:Y-1] - dx_ZUx[1:X-1,1:Y-1]

    fn_p[1:X-1,1:Y-1]  = CIP(coeff31[0],  coeff31[1],  coeff31[2],  coeff31[3], \
                               f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                               g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]  = CIP(coeff31[4],  coeff31[5],  coeff31[6],  coeff31[7], \
                               f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                               g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]  = CIP(coeff31[8],  coeff31[9],  coeff31[10], coeff31[11], \
                               f_m[2:X,1:Y-1],   f_m[1:X-1,1:Y-1], \
                               g_m[2:X,1:Y-1],   g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]  = CIP(coeff31[12], coeff31[13], coeff31[14], coeff31[15], \
                               f_m[2:X,1:Y-1],   f_m[1:X-1,1:Y-1], \
                               g_m[2:X,1:Y-1],   g_m[1:X-1,1:Y-1])

    P[1:X-1,1:Y-1]     = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) 
    ZUx[1:X-1,1:Y-1]    = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) 
    dx_P[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) 
    dx_ZUx[1:X-1,1:Y-1] = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) 

    f_p[1:X-1,1:Y-1]   = dx_P[1:X-1,1:Y-1] \
                          + (ZUy[2:X,1:Y-1] - ZUy[0:X-2,1:Y-1]) / (2. * dx)
    f_m[1:X-1,1:Y-1]   = dx_P[1:X-1,1:Y-1] \
                          - (ZUy[2:X,1:Y-1] - ZUy[0:X-2,1:Y-1]) / (2. * dx)

    fn_p[1:X-1,1:Y-1]  = LINEAR(coeff10[4], coeff10[5], \
                                 f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]  = LINEAR(coeff10[6], coeff10[7], \
                                 f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1])

    dx_P[1:X-1,1:Y-1]  = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) 

    f_p[1:X-1,1:Y-1]   =    P[1:X-1,1:Y-1] + ZUy[1:X-1,1:Y-1]
    f_m[1:X-1,1:Y-1]   =    P[1:X-1,1:Y-1] -   ZUy[1:X-1,1:Y-1]
    g_p[1:X-1,1:Y-1]   = dy_P[1:X-1,1:Y-1] + dy_ZUy[1:X-1,1:Y-1]
    g_m[1:X-1,1:Y-1]   = dy_P[1:X-1,1:Y-1] - dy_ZUy[1:X-1,1:Y-1]

    fn_p[1:X-1,1:Y-1]  = CIP(coeff31[16], coeff31[17], coeff31[18], coeff31[19], \
                               f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1], \
                               g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]  = CIP(coeff31[20], coeff31[21], coeff31[22], coeff31[23], \
                               f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1], \
                               g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]  = CIP(coeff31[24], coeff31[25], coeff31[26], coeff31[27], \
                               f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                               g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]  = CIP(coeff31[28], coeff31[29], coeff31[30], coeff31[31], \
                               f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                               g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])

    P[1:X-1,1:Y-1]     = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) 
    ZUy[1:X-1,1:Y-1]    = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) 
    dy_P[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) 
    dy_ZUy[1:X-1,1:Y-1] = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) 

end = clock()

print ("Processing Time : " + str((end - start)) + " [sec]")

