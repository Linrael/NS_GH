import numpy as np
from matplotlib import pyplot as plt

SMALL_SIZE=30
MEDIUM_SIZE =SMALL_SIZE+4
BIGGER_SIZE = SMALL_SIZE+10
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


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

startx = 0
endx = 6
starty = 0
endy = 0.5

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

plot_number=0

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

    indx=5
    indy=5

    maskedP = np.ma.array(P, mask=fluid_flag)
    maskedP = maskedP[1:-1, 1:-1]
    maskedP = maskedP[::-1, :]+6.2
    #maskedP=np.exp(maskedP)

    if plot_number  == 8:
        fig,ax= plt.subplots(figsize=(30,20))

        color=np.sqrt(np.hypot(U,V))
        stream=ax.streamplot(X, Y, U, V, color=color, linewidth=4,cmap='afmhot',arrowstyle='-|>',arrowsize=2)
        ax.set_title(f'Backwards-Facing Step at 11 seconds')
        ax.imshow(konvertieren,alpha=0.8, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax],
                  interpolation='nearest')
        im=ax.imshow(maskedP, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax],alpha=0.99,norm="linear",vmin=-0.3,vmax=0.35)
        ax.set_xlim(0, 3.5)
        ax.set_ylim(0, 0.5)
        # ax.set_xlabel("x")
        # ax.set_ylabel("y")
        #
        fig.colorbar(im,ax=ax,label="Relative pressure P",extend='both', shrink = 0.37, fraction = 0.1,orientation='horizontal',location='bottom',anchor=(0.85,0.),pad=-0.12)
        fig.colorbar(stream.lines,ax=ax,label="$\Vert u \Vert $",location="bottom",shrink=0.37,orientation='horizontal',anchor=(0.15,0.15),pad=0.05)
        #fig.subplots_adjust(left=0.05, right=0.95, bottom=0.05, top=0.95, hspace=0.1,wspace=0.1)
        ax.set_aspect(3)
        plt.savefig('bwstep11.pdf',dpi=300.)
        plt.show()
    plot_number = plot_number + 1




# while True:
#     line = f.readline()
#     if not line:
#         break
#     timestep = int(line)
#     U = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
#     V = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
#     P = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
#     # choose choatic part:
#     U = U[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
#     V = V[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
#     maskedP = np.ma.array(P, mask=fluid_flag)
#     maskedP = maskedP[1:-1, :]
#     maskedP = maskedP[int(starty / ylength * jmax) + 1: int(endy / ylength * jmax) + 1, int(startx / xlength * imax) + 1: int(endx / xlength * imax) + 1]
#     maskedP = maskedP[::-1, :]
#     #
#     print(V.shape)
#
#     fig, ax = plt.subplots(1)
#     ax.streamplot(Xsmall, Ysmall, U, V, color=U, linewidth=0.5, cmap='autumn')
#     ax.set_title(f'U V Stream Plot at Timestep {timestep} for time# {timestep * delt}')
#     ax.imshow(konvertieren, extent=[startx * xlength, endx * xlength, starty * ylength, endy * ylength], interpolation='nearest')
#     ax.imshow(maskedP, extent=[startx, endx, starty, endy])
#     plt.show()

    # 11.036 seconds is last data for bs
# xlength 6.
# ylength 0.5
# imax 1800.
# jmax 450.
# delx 0.01
# dely 0.01
# t_end 20.
# delt 0.02
# tau 0.5
# N 0.
# itermax 7500.
# eps 0.001
# omega 1.7
# gamma 0.9
# Re 500.
# Pr 1.0
# beta 0.
# GX 0.
# GY 0.
# UI 0.
# VI 0.
# PI 0.
