# 3 dimension cip acoustic simulation

from time  import *
from pylab import *
from numpy import *
from scipy import *
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

X = 101
Y = 101
NT = 601

f_p    = zeros((X,Y), "float64")
f_m    = zeros((X,Y), "float64")
g_p    = zeros((X,Y), "float64")
g_m    = zeros((X,Y), "float64")
fn_p   = zeros((X,Y), "float64")
fn_m   = zeros((X,Y), "float64")
gn_p   = zeros((X,Y), "float64")
gn_m   = zeros((X,Y), "float64")

P      = zeros((X,Y), "float64")
dx_P   = zeros((X,Y), "float64")
dy_P   = zeros((X,Y), "float64")
dxy_P  = zeros((X,Y), "float64")

Ux     = zeros((X,Y), "float64")
dx_Ux  = zeros((X,Y), "float64")
dy_Ux  = zeros((X,Y), "float64")
dxy_Ux = zeros((X,Y), "float64")

Uy     = zeros((X,Y), "float64")
dx_Uy  = zeros((X,Y), "float64")
dy_Uy  = zeros((X,Y), "float64")
dxy_Uy = zeros((X,Y), "float64")

xc= (X-1) / 2
yc= (Y-1) / 2

dx = 5.e-2
dy = 5.e-2
dt = 5.e-5

Ro = 1.21
bm = 1.4235529e5
c0 = sqrt(bm / Ro)
Z0 = sqrt(bm * Ro)
sigma = 0.2

mic = zeros((NT,4), "float64")

coeff  = zeros(32, "float64")

Ua = c0
xi =-Ua * dt
C  = c0 * dt / dx
C2 = C  * C
C3 = C2 * C

coeff[0]  = (-2. * C3 + 3. * C2)
coeff[1]  = (2. * C3 - 3. * C2 + 1.)
coeff[2]  = xi * (C2 - C)
coeff[3]  = xi * (C2 - 2. * C + 1.)
coeff[4]  = 6. * (-C3 + C2) / xi
coeff[5]  = 6. * (C3 - C2) / xi
coeff[6]  = (3. * C2 - 2. * C)
coeff[7]  = (3. * C2 - 4. * C + 1.)

Ua =-c0
xi =-Ua * dt
C  = c0 * dt / dx
C2 = C  * C
C3 = C2 * C

coeff[8]  = (-2. * C3 + 3. * C2)
coeff[9]  = (2. * C3 - 3. * C2 + 1.)
coeff[10] = xi * (C2 - C)
coeff[11] = xi * (C2 - 2. * C + 1.)
coeff[12] = 6. * (-C3 + C2) / xi
coeff[13] = 6. * (C3 - C2) / xi
coeff[14] = (3. * C2 - 2. * C)
coeff[15] = (3. * C2 - 4. * C + 1.)

Ua = c0
yi =-Ua * dt
C  = c0 * dt / dy
C2 = C  * C
C3 = C2 * C

coeff[16] = (-2. * C3 + 3. * C2)
coeff[17] = (2. * C3 - 3. * C2 + 1.)
coeff[18] = yi * (C2 - C)
coeff[19] = yi * (C2 - 2. * C + 1.)
coeff[20] = 6. * (-C3 + C2) / yi
coeff[21] = 6. * (C3 - C2) / yi
coeff[22] = (3. * C2 - 2. * C)
coeff[23] = (3. * C2 - 4. * C + 1.)

Ua =-c0
yi =-Ua * dt
C  = c0 * dt / dy
C2 = C  * C
C3 = C2 * C

coeff[24] = (-2. * C3 + 3. * C2)
coeff[25] = (2. * C3 - 3. * C2 + 1.)
coeff[26] = yi * (C2 - C)
coeff[27] = yi * (C2 - 2. * C + 1.)
coeff[28] = 6. * (-C3 + C2) / yi
coeff[29] = 6. * (C3 - C2) / yi
coeff[30] = (3. * C2 - 2. * C)
coeff[31] = (3. * C2 - 4. * C + 1.)

for i in range(1, X-1):
    x = dx * i
    for j in range(1, Y-1):
        y = dy * j
        TX = x - xc * dx
        TY = y - yc * dy
        P[i][j]     = \
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma**2))
        dx_P[i][j]  = -TX * \
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma**2)) / sigma**2
        dy_P[i][j]  = -TY * \
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma**2)) / sigma**2
        dxy_P[i][j] = TX * TY * \
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma**2)) / sigma**4

def CIP(coeff0, coeff1, coeff2, coeff3, f0, f1, g0, g1):
    return    coeff0 * f0 \
            + coeff1 * f1 \
            + coeff2 * g0 \
            + coeff3 * g1

start = clock()

for t in range(NT):
    if (t % 50) == 0:
        figure(figsize=(7.5,6))
        pcolor(P)
        colorbar()
        xlim(1,X-1)
        ylim(1,Y-1)
        clim(-0.25,0.25)
        xlabel("Y [sample]")
        ylabel("X [sample]")
        show()
        
    RX=0.5
    RY=0.5
    
    f_p[1:X,  1:Y-1]    =    P[1:X,  1:Y-1] + (Z0 *    Ux[1:X,  1:Y-1])
    g_p[1:X,  1:Y-1]    = dx_P[1:X,  1:Y-1] + (Z0 * dx_Ux[1:X,  1:Y-1])
    f_m[0:X-1,1:Y-1]    =    P[0:X-1,1:Y-1] - (Z0 *    Ux[0:X-1,1:Y-1])
    g_m[0:X-1,1:Y-1]    = dx_P[0:X-1,1:Y-1] - (Z0 * dx_Ux[0:X-1,1:Y-1])

    fn_p[2:X  ,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                f_p[1:X-1,1:Y-1], f_p[2:X  ,1:Y-1], \
                                g_p[1:X-1,1:Y-1], g_p[2:X  ,1:Y-1])
    gn_p[2:X  ,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                f_p[1:X-1,1:Y-1], f_p[2:X  ,1:Y-1], \
                                g_p[1:X-1,1:Y-1], g_p[2:X  ,1:Y-1])
    fn_m[0:X-2,1:Y-1]   = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                f_m[1:X-1,1:Y-1], f_m[0:X-2,1:Y-1], \
                                g_m[1:X-1,1:Y-1], g_m[0:X-2,1:Y-1])
    gn_m[0:X-2,1:Y-1]   = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                f_m[1:X-1,1:Y-1], f_m[0:X-2,1:Y-1], \
                                g_m[1:X-1,1:Y-1], g_m[0:X-2,1:Y-1])

    fn_p[1:2  ,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                RX*f_m[0:1,1:Y-1], f_p[1:2  ,1:Y-1], \
                                -RX*g_m[0:1,1:Y-1], g_p[1:2  ,1:Y-1])
    gn_p[1:2  ,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                RX*f_m[0:1,1:Y-1], f_p[1:2  ,1:Y-1], \
                                -RX*g_m[0:1,1:Y-1], g_p[1:2  ,1:Y-1])
    fn_m[X-2:X-1,1:Y-1] = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                RX*f_p[X-1:X,1:Y-1], f_m[X-2:X-1,  1:Y-1], \
                                -RX*g_p[X-1:X,1:Y-1], g_m[X-2:X-1,1:Y-1])
    gn_m[X-2:X-1,1:Y-1] = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                RX*f_p[X-1:X,1:Y-1], f_m[X-2:X-1,  1:Y-1], \
                                -RX*g_p[X-1:X,1:Y-1], g_m[X-2:X-1,1:Y-1])

    P[:  ,1:Y-1]      = (fn_p[:  ,1:Y-1] + fn_m[:  ,1:Y-1]) /  2.
    Ux[:  ,1:Y-1]     = (fn_p[:  ,1:Y-1] - fn_m[:  ,1:Y-1]) / (2. * Z0)
    dx_P[:  ,1:Y-1]   = (gn_p[:  ,1:Y-1] + gn_m[:  ,1:Y-1]) /  2.
    dx_Ux[:  ,1:Y-1]  = (gn_p[:  ,1:Y-1] - gn_m[:  ,1:Y-1]) / (2. * Z0)


    f_p[1:X  ,1:Y-1]    =  dy_P[1:X  ,1:Y-1] + (Z0 *  dy_Ux[1:X  ,1:Y-1])
    g_p[1:X  ,1:Y-1]    = dxy_P[1:X  ,1:Y-1] + (Z0 * dxy_Ux[1:X  ,1:Y-1])
    f_m[0:X-1,1:Y-1]    =  dy_P[0:X-1,1:Y-1] - (Z0 *  dy_Ux[0:X-1,1:Y-1])
    g_m[0:X-1,1:Y-1]    = dxy_P[0:X-1,1:Y-1] - (Z0 * dxy_Ux[0:X-1,1:Y-1])

    fn_p[2:X  ,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                f_p[1:X-1,1:Y-1], f_p[2:X  ,1:Y-1], \
                                g_p[1:X-1,1:Y-1], g_p[2:X  ,1:Y-1])
    gn_p[2:X  ,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                f_p[1:X-1,1:Y-1], f_p[2:X  ,1:Y-1], \
                                g_p[1:X-1,1:Y-1], g_p[2:X  ,1:Y-1])
    fn_m[0:X-2,1:Y-1]   = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                f_m[1:X-1,1:Y-1], f_m[0:X-2,1:Y-1], \
                                g_m[1:X-1,1:Y-1], g_m[0:X-2,1:Y-1])
    gn_m[0:X-2,1:Y-1]   = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                f_m[1:X-1,1:Y-1], f_m[0:X-2,1:Y-1], \
                                g_m[1:X-1,1:Y-1], g_m[0:X-2,1:Y-1])

    fn_p[1:2  ,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                RX*f_m[0:1,1:Y-1], f_p[1:2  ,1:Y-1], \
                                -RX*g_m[0:1,1:Y-1], g_p[1:2  ,1:Y-1])
    gn_p[1:2  ,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                RX*f_m[0:1,1:Y-1], f_p[1:2  ,1:Y-1], \
                                -RX*g_m[0:1,1:Y-1], g_p[1:2  ,1:Y-1])
    fn_m[X-2:X-1,1:Y-1] = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                RX*f_p[X-1:X,1:Y-1], f_m[X-2:X-1,  1:Y-1], \
                                -RX*g_p[X-1:X,1:Y-1], g_m[X-2:X-1,1:Y-1])
    gn_m[X-2:X-1,1:Y-1] = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                RX*f_p[X-1:X,1:Y-1], f_m[X-2:X-1,  1:Y-1], \
                                -RX*g_p[X-1:X,1:Y-1], g_m[X-2:X-1,1:Y-1])

    dy_P[:  ,1:Y-1]   = (fn_p[:  ,1:Y-1] + fn_m[:  ,1:Y-1]) /  2.
    dy_Ux[:  ,1:Y-1]  = (fn_p[:  ,1:Y-1] - fn_m[:  ,1:Y-1]) / (2. * Z0)
    dxy_P[:  ,1:Y-1]  = (gn_p[:  ,1:Y-1] + gn_m[:  ,1:Y-1]) /  2.
    dxy_Ux[:  ,1:Y-1] = (gn_p[:  ,1:Y-1] - gn_m[:  ,1:Y-1]) / (2. * Z0)


    f_p[1:X-1,1:Y]    =    P[1:X-1,1:Y] + (Z0 *    Uy[1:X-1,1:Y])
    g_p[1:X-1,1:Y]    = dy_P[1:X-1,1:Y] + (Z0 * dy_Uy[1:X-1,1:Y])
    f_m[1:X-1,0:Y-1]    =    P[1:X-1,0:Y-1] - (Z0 *    Uy[1:X-1,0:Y-1])
    g_m[1:X-1,0:Y-1]    = dy_P[1:X-1,0:Y-1] - (Z0 * dy_Uy[1:X-1,0:Y-1])

    fn_p[1:X-1,2:Y]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                              f_p[1:X-1,1:Y-1], f_p[1:X-1,2:Y], \
                              g_p[1:X-1,1:Y-1], g_p[1:X-1,2:Y])
    gn_p[1:X-1,2:Y]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                              f_p[1:X-1,1:Y-1], f_p[1:X-1,2:Y], \
                              g_p[1:X-1,1:Y-1], g_p[1:X-1,2:Y])
    fn_m[1:X-1,0:Y-2] = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                              f_m[1:X-1,1:Y-1],   f_m[1:X-1,0:Y-2], \
                              g_m[1:X-1,1:Y-1],   g_m[1:X-1,0:Y-2])
    gn_m[1:X-1,0:Y-2] = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                              f_m[1:X-1,1:Y-1],   f_m[1:X-1,0:Y-2], \
                              g_m[1:X-1,1:Y-1],   g_m[1:X-1,0:Y-2])

    fn_p[1:X-1,1:2]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                              RY*f_m[1:X-1,0:1], f_p[1:X-1,1:2], \
                              -RY*g_m[1:X-1,0:1], g_p[1:X-1,1:2])
    gn_p[1:X-1,1:2]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                              RY*f_m[1:X-1,0:1], f_p[1:X-1,1:2], \
                              -RY*g_m[1:X-1,0:1], g_p[1:X-1,1:2])
    
    fn_m[1:X-1,Y-2:Y-1]   = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                                  RY*f_p[1:X-1,Y-1:Y],   f_m[1:X-1,Y-2:Y-1], \
                                  -RY*g_p[1:X-1,Y-1:Y],   g_m[1:X-1,Y-2:Y-1])
    gn_m[1:X-1,Y-2:Y-1]   = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                                  RY*f_p[1:X-1,Y-1:Y],   f_m[1:X-1,Y-2:Y-1], \
                                  -RY*g_p[1:X-1,Y-1:Y],   g_m[1:X-1,Y-2:Y-1])

    P[1:X-1,0:Y]      = (fn_p[1:X-1,0:Y] + fn_m[1:X-1,0:Y]) /  2.
    Uy[1:X-1,0:Y]     = (fn_p[1:X-1,0:Y] - fn_m[1:X-1,0:Y]) / (2. * Z0)
    dy_P[1:X-1,0:Y]   = (gn_p[1:X-1,0:Y] + gn_m[1:X-1,0:Y]) /  2.
    dy_Uy[1:X-1,0:Y]  = (gn_p[1:X-1,0:Y] - gn_m[1:X-1,0:Y]) / (2. * Z0)


    f_p[1:X-1,1:Y]    =  dx_P[1:X-1,1:Y] + (Z0 *  dx_Uy[1:X-1,1:Y])
    g_p[1:X-1,0:Y]    = dxy_P[1:X-1,0:Y] + (Z0 * dxy_Uy[1:X-1,0:Y])
    f_m[1:X-1,0:Y-1]    =  dx_P[1:X-1,0:Y-1] - (Z0 *  dx_Uy[1:X-1,0:Y-1])
    g_m[1:X-1,0:Y-1]    = dxy_P[1:X-1,0:Y-1] - (Z0 * dxy_Uy[1:X-1,0:Y-1])

    fn_p[1:X-1,2:Y]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                              f_p[1:X-1,1:Y-1], f_p[1:X-1,2:Y], \
                              g_p[1:X-1,1:Y-1], g_p[1:X-1,2:Y])
    gn_p[1:X-1,2:Y]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                              f_p[1:X-1,1:Y-1], f_p[1:X-1,2:Y], \
                              g_p[1:X-1,1:Y-1], g_p[1:X-1,2:Y])
    fn_m[1:X-1,0:Y-2] = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                              f_m[1:X-1,1:Y-1],   f_m[1:X-1,0:Y-2], \
                              g_m[1:X-1,1:Y-1],   g_m[1:X-1,0:Y-2])
    gn_m[1:X-1,0:Y-2] = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                              f_m[1:X-1,1:Y-1],   f_m[1:X-1,0:Y-2], \
                              g_m[1:X-1,1:Y-1],   g_m[1:X-1,0:Y-2])

    fn_p[1:X-1,1:2]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                              RY*f_m[1:X-1,0:1], f_p[1:X-1,1:2], \
                              -RY*g_m[1:X-1,0:1], g_p[1:X-1,1:2])
    gn_p[1:X-1,1:2]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                              RY*f_m[1:X-1,0:1], f_p[1:X-1,1:2], \
                              -RY*g_m[1:X-1,0:1], g_p[1:X-1,1:2])
    fn_m[1:X-1,Y-2:Y-1]   = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                                  RY*f_p[1:X-1,Y-1:Y],   f_m[1:X-1,Y-2:Y-1], \
                                  -RY*g_p[1:X-1,Y-1:Y],   g_m[1:X-1,Y-2:Y-1])
    gn_m[1:X-1,Y-2:Y-1]   = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                                  RY*f_p[1:X-1,Y-1:Y],   f_m[1:X-1,Y-2:Y-1], \
                                  -RY*g_p[1:X-1,Y-1:Y],   g_m[1:X-1,Y-2:Y-1])

    dx_P[1:X-1,0:Y]   = (fn_p[1:X-1,0:Y] + fn_m[1:X-1,0:Y]) /  2.
    dx_Uy[1:X-1,0:Y]  = (fn_p[1:X-1,0:Y] - fn_m[1:X-1,0:Y]) / (2. * Z0)
    dxy_P[1:X-1,0:Y]  = (gn_p[1:X-1,0:Y] + gn_m[1:X-1,0:Y]) /  2.
    dxy_Uy[1:X-1,0:Y] = (gn_p[1:X-1,0:Y] - gn_m[1:X-1,0:Y]) / (2. * Z0)


end = clock()

print "Processing Time : " + str((end - start)) + " [sec]"

