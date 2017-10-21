# 1 dimension cip acoustic simulation
from time  import *
from pylab import *
from numpy import *
from scipy import *
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d

X = 201
Y = 201
NT = 501

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

mic = zeros((NT,3), "float64")

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

f = open('mic.txt', 'w') 

def CIP(coeff0, coeff1, coeff2, coeff3, f0, f1, g0, g1):
    return    coeff0 * f0 \
            + coeff1 * f1 \
            + coeff2 * g0 \
            + coeff3 * g1

start = clock()

for t in range(NT):

    mic[t,0] = P[xc+20,yc]
    mic[t,1] = P[xc+30,yc]
    mic[t,2] = P[xc+40,yc]

    f.write(str(t*dt)+"\t")
    f.write(str(mic[t,0])+"\t")
    f.write(str(mic[t,1])+"\t")
    f.write(str(mic[t,2])+"\n")    

    f_p[1:X-1,1:Y-1]    =    P[1:X-1,1:Y-1] + (Z0 *    Ux[1:X-1,1:Y-1])
    f_m[1:X-1,1:Y-1]    =    P[1:X-1,1:Y-1] - (Z0 *    Ux[1:X-1,1:Y-1])
    g_p[1:X-1,1:Y-1]    = dx_P[1:X-1,1:Y-1] + (Z0 * dx_Ux[1:X-1,1:Y-1])
    g_m[1:X-1,1:Y-1]    = dx_P[1:X-1,1:Y-1] - (Z0 * dx_Ux[1:X-1,1:Y-1])

    fn_p[1:X-1,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                                g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                                g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]   = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                f_m[2:X,1:Y-1], f_m[1:X-1,1:Y-1], \
                                g_m[2:X,1:Y-1], g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]   = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                f_m[2:X,1:Y-1], f_m[1:X-1,1:Y-1], \
                                g_m[2:X,1:Y-1], g_m[1:X-1,1:Y-1])

    P[1:X-1,1:Y-1]      = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) /  2.
    Ux[1:X-1,1:Y-1]     = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) / (2. * Z0)
    dx_P[1:X-1,1:Y-1]   = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) /  2.
    dx_Ux[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) / (2. * Z0)

    f_p[1:X-1,1:Y-1]    =  dy_P[1:X-1,1:Y-1] + (Z0 *  dy_Ux[1:X-1,1:Y-1])
    f_m[1:X-1,1:Y-1]    =  dy_P[1:X-1,1:Y-1] - (Z0 *  dy_Ux[1:X-1,1:Y-1])
    g_p[1:X-1,1:Y-1]    = dxy_P[1:X-1,1:Y-1] + (Z0 * dxy_Ux[1:X-1,1:Y-1])
    g_m[1:X-1,1:Y-1]    = dxy_P[1:X-1,1:Y-1] - (Z0 * dxy_Ux[1:X-1,1:Y-1])

    fn_p[1:X-1,1:Y-1]   = CIP(coeff[0],  coeff[1],  coeff[2],  coeff[3], \
                                f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                                g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]   = CIP(coeff[4],  coeff[5],  coeff[6],  coeff[7], \
                                f_p[0:X-2,1:Y-1], f_p[1:X-1,1:Y-1], \
                                g_p[0:X-2,1:Y-1], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]   = CIP(coeff[8],  coeff[9],  coeff[10], coeff[11], \
                                f_m[2:X,1:Y-1],   f_m[1:X-1,1:Y-1], \
                                g_m[2:X,1:Y-1],   g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]   = CIP(coeff[12], coeff[13], coeff[14], coeff[15], \
                                f_m[2:X,1:Y-1],   f_m[1:X-1,1:Y-1], \
                                g_m[2:X,1:Y-1],   g_m[1:X-1,1:Y-1])

    dy_P[1:X-1,1:Y-1]   = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) /  2.
    dy_Ux[1:X-1,1:Y-1]  = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) / (2. * Z0)
    dxy_P[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) /  2.
    dxy_Ux[1:X-1,1:Y-1] = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) / (2. * Z0)

    f_p[1:X-1,1:Y-1]    =    P[1:X-1,1:Y-1] + (Z0 *    Uy[1:X-1,1:Y-1])
    f_m[1:X-1,1:Y-1]    =    P[1:X-1,1:Y-1] - (Z0 *    Uy[1:X-1,1:Y-1])
    g_p[1:X-1,1:Y-1]    = dy_P[1:X-1,1:Y-1] + (Z0 * dy_Uy[1:X-1,1:Y-1])
    g_m[1:X-1,1:Y-1]    = dy_P[1:X-1,1:Y-1] - (Z0 * dy_Uy[1:X-1,1:Y-1])

    fn_p[1:X-1,1:Y-1]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                                f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1],
                                g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                                f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1], \
                                g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]   = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                                f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                                g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]   = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                                f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                                g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])

    P[1:X-1,1:Y-1]      = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) /  2.
    Uy[1:X-1,1:Y-1]     = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) / (2. * Z0)
    dy_P[1:X-1,1:Y-1]   = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) /  2.
    dy_Uy[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) / (2. * Z0)

    f_p[1:X-1,1:Y-1]    =  dx_P[1:X-1,1:Y-1] + (Z0 *  dx_Uy[1:X-1,1:Y-1])
    f_m[1:X-1,1:Y-1]    =  dx_P[1:X-1,1:Y-1] - (Z0 *  dx_Uy[1:X-1,1:Y-1])
    g_p[1:X-1,1:Y-1]    = dxy_P[1:X-1,1:Y-1] + (Z0 * dxy_Uy[1:X-1,1:Y-1])
    g_m[1:X-1,1:Y-1]    = dxy_P[1:X-1,1:Y-1] - (Z0 * dxy_Uy[1:X-1,1:Y-1])

    fn_p[1:X-1,1:Y-1]   = CIP(coeff[16], coeff[17], coeff[18], coeff[19], \
                                f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1], \
                                g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    gn_p[1:X-1,1:Y-1]   = CIP(coeff[20], coeff[21], coeff[22], coeff[23], \
                                f_p[1:X-1,0:Y-2], f_p[1:X-1,1:Y-1], \
                                g_p[1:X-1,0:Y-2], g_p[1:X-1,1:Y-1])
    fn_m[1:X-1,1:Y-1]   = CIP(coeff[24], coeff[25], coeff[26], coeff[27], \
                                f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                                g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])
    gn_m[1:X-1,1:Y-1]   = CIP(coeff[28], coeff[29], coeff[30], coeff[31], \
                                f_m[1:X-1,2:Y],   f_m[1:X-1,1:Y-1], \
                                g_m[1:X-1,2:Y],   g_m[1:X-1,1:Y-1])

    dx_P[1:X-1,1:Y-1]   = (fn_p[1:X-1,1:Y-1] + fn_m[1:X-1,1:Y-1]) /  2.
    dx_Uy[1:X-1,1:Y-1]  = (fn_p[1:X-1,1:Y-1] - fn_m[1:X-1,1:Y-1]) / (2. * Z0)
    dxy_P[1:X-1,1:Y-1]  = (gn_p[1:X-1,1:Y-1] + gn_m[1:X-1,1:Y-1]) /  2.
    dxy_Uy[1:X-1,1:Y-1] = (gn_p[1:X-1,1:Y-1] - gn_m[1:X-1,1:Y-1]) / (2. * Z0)


figure(1)
plot(mic)
xlabel("Sampling point")
ylabel("Relative sound pressure")
show()

f.close()
end = clock()

print "Processing Time : " + str((end - start)) + " [sec]"

