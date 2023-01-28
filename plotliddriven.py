import numpy as np
from matplotlib import pyplot as plt

SMALL_SIZE=36
MEDIUM_SIZE =SMALL_SIZE+4
BIGGER_SIZE = SMALL_SIZE+8
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


f = open('finaldata/liddrivencavityplot.txt', 'r')

# read parameters

imax = int(f.readline())
jmax = int(f.readline())
xlength = float(f.readline())
ylength = float(f.readline())
delt = float(f.readline())
delt=0.001
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

malfig, malax = plt.subplots()
malax.imshow(konvertieren, extent=[0, 3, 0, 0.5], interpolation='nearest')
malax.set_title("Geometry")
plt.show()

fluid_flag = np.logical_not(fluid_flag)

X = np.arange(0, xlength + 0 * 1.1 * xlength / imax, xlength / imax)
Y = np.arange(0, ylength + 0 * 1.1 * ylength / jmax, ylength / jmax)
X, Y = np.meshgrid(X, Y)

plot_number = 0

while True:
    line = f.readline()
    if not line:
        break
    timestep = int(line)
    U = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    V = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    P = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)

    U = U[1:-1, 1:-1]
    V = V[1:-1, 1:-1]

    maskedP = np.ma.array(P, mask=fluid_flag)
    maskedP = maskedP[1:-1, 1:-1]
    maskedP = maskedP[::-1, :]
    #maskedP=np.exp(maskedP)

    if plot_number % 1 == 0:
        fig,ax= plt.subplots(figsize=(20,18))
        stream=ax.streamplot(X, Y, U, V, color=U, linewidth=3, cmap='autumn')
        ax.set_title(f'Lid-Driven Cavity at {timestep * delt} seconds')
        ax.imshow(konvertieren, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax],
                  interpolation='nearest')
        im=ax.imshow(maskedP, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax],alpha=0.99,norm="linear",vmax=0.15,vmin=-0.07)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        fig.colorbar(im,ax=ax,label="Relative pressure P",shrink=0.5)
        fig.colorbar(stream.lines,ax=ax,label="Velocity field",location="left",shrink=0.5,pad=0.13)
        plt.show()

    plot_number = plot_number + 1