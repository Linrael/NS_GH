import numpy as np
from matplotlib import pyplot as plt

f = open('finaldata/smallbackwardstep.txt', 'r')

# read parameters


imax = int(f.readline())
jmax = int(f.readline())
xlength = float(f.readline())
ylength = float(f.readline())
delt = float(f.readline())
print(imax, jmax, xlength, ylength, delt, sep=" ")

# initial fields

U = f.readline().split('/')[:-1]
U = np.array(U)
U = 1 * U.astype(np.double).reshape(imax + 2, jmax + 2)
U = np.swapaxes(U, 0, 1)

V = f.readline().split('/')[:-1]
V = np.array(V)
V = 1 * V.astype(np.double).reshape(imax + 2, jmax + 2)
V = np.swapaxes(V, 0, 1)

P = f.readline().split('/')[:-1]
P = np.array(P)
P = P.astype(np.double).reshape(imax + 2, jmax + 2)
P = np.swapaxes(P, 0, 1)

# geometry

FLAG = f.readline().split('/')[:-1]
FLAG = np.array(FLAG).astype(np.double)
fluid_flag = [1. if x == 16 else 0. for x in FLAG]
fluid_flag = np.array(fluid_flag).reshape(imax + 2, jmax + 2).astype(np.double)
fluid_flag = np.swapaxes(fluid_flag, 0, 1)
malen = fluid_flag[1:-1, 1:-1]  # remove for ghostpoints
konvertieren = np.zeros((jmax, imax, 3))  # +2 for ghost points
konvertieren[malen > 0.5] = [1, 1, 1]
konvertieren[malen < 0.5] = [0, 0, 0]
konvertieren = konvertieren[::-1, :]

# malfig, malax = plt.subplots()
# malax.imshow(konvertieren, extent=[0, 3, 0, 0.5], interpolation='nearest')
# malax.set_title("Geometry")
# plt.show()

fluid_flag = np.logical_not(fluid_flag)

startx = 1.75
endx = 1.8
starty = 0.1
endy = 0.12

X = np.arange(0, xlength + 0 * 1.1 * xlength / imax, xlength / imax)
Y = np.arange(0, ylength + 0 * 1.1 * ylength / jmax, ylength / jmax)
X, Y = np.meshgrid(X, Y)
# Xsmall2 = X[int(startx * imax):int(endx * imax), int(starty * jmax):int(endy * jmax)]
# Ysmall = Y[int(startx * imax):int(endx * imax), int(starty * jmax):int(endy * jmax)]
# print(Xsmall2.shape)

Xsmall = np.arange(startx, endx, xlength/imax)
Ysmall = np.arange(starty, endy, ylength/jmax)
Xsmall, Ysmall = np.meshgrid(Xsmall, Ysmall)
print(Ysmall.shape)
# now read additional data
while True:
    line = f.readline()
    if not line:
        break
    timestep = int(line)
    U = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    V = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    P = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    # choose choatic part:
    U = U[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
    V = V[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
    # maskedP = np.ma.array(P, mask=fluid_flag)
    # maskedP = maskedP[1:-1, :]
    # maskedP = maskedP[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
    # maskedP = maskedP[::-1, :]
    #
    print(V.shape)

    fig, ax = plt.subplots(1)
    ax.streamplot(Xsmall, Ysmall, U, V, color=U, linewidth=0.5, cmap='autumn')
    ax.set_title(f'U V Stream Plot at Timestep {timestep} for time# {timestep * delt}')
    # ax.imshow(konvertieren, extent=[startx * xlength, endx * xlength, starty * ylength, endy * ylength], interpolation='nearest')
    # ax.imshow(maskedP, extent=[startx, endx, starty, endy])
    plt.show()

    # 11.036 seconds is last data for bs
