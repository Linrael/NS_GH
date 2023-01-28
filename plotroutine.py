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

X = np.arange(0, xlength + 0 * 1.1 * xlength / imax, xlength / imax)
Y = np.arange(0, ylength + 0 * 1.1 * ylength / jmax, ylength / jmax)
X, Y = np.meshgrid(X, Y)
# X = np.ma.array(X, mask=fluid_flag[1:-1, 1:-1])
# Y = np.ma.array(Y, mask=fluid_flag[1:-1, 1:-1])

# fig, ax = plt.subplots()
# ax.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
# plt.show()
#
# fig, ax = plt.subplots()
# ax.quiver(X, Y, U, V)
# ax.set_title('U V Plot')
# plt.show()

# now read additional data
indx = 5  # jeder 4. datenpunkte wird geplottet
indy = 5
plot_number = 0

while True:
    line = f.readline()
    if not line:
        break
    timestep = int(line)
    U = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    V = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    P = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)

    # print(timestep, np.sum(U[:, 1] + U[:,2]), np.sum(U[:, -2] + U[:,-3]))

    maskedU = np.ma.array(U, mask=fluid_flag)
    maskedU = maskedU[1:-1, 1:-1]
    maskedV = np.ma.array(V, mask=fluid_flag)
    maskedV = maskedV[1:-1, 1:-1]
    U = U[1:-1, 1:-1]
    V = V[1:-1, 1:-1]
    maskedP = np.ma.array(P, mask=fluid_flag)
    maskedP = maskedP[1:-1, 1:-1]
    maskedP = maskedP[::-1, :]
    if plot_number % 3 == 0:
        fig, ax = plt.subplots()
        ax.streamplot(X, Y, U, V, color=U, linewidth=.5, cmap='autumn')
        ax.set_title(f'U V Stream Plot at Timestep {timestep} for time# {timestep * delt}')
        ax.imshow(konvertieren, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax], interpolation='nearest')
        ax.imshow(maskedP, extent=[0, xlength - xlength / imax, 0, ylength - ylength / jmax])
        ax.spines['right'].set_color('none')
        ax.spines['left'].set_color('none')
        plt.show()
    break
        # fig, ax = plt.subplots()
        # ax.quiver(X[::indx, ::indy], Y[::indx, ::indy], U[::indx, ::indy], V[::indx, ::indy])
        # ax.set_title(f'U V Quiver Plot at Timestep {timestep} for time {timestep * delt}')
        # plt.show()
    plot_number = plot_number + 1



# to save all data use this:
# U=[]
# V=[]
# P=[]
# timesteps=[]
#
#

#
# for i in range(len(timesteps)):
#     ax[2*i].streamplot(X, Y, U[0], V[0], color=U[0], linewidth=2, cmap='autumn')
#     ax[2*i].set_title(f'U V Stream Plot at {timesteps[i]}')
#     ax[2*i+1].quiver(X, Y, U[0], V[0])
#     ax[2*i+1].set_title(f'U V Quiver Plot at {timesteps[i]}')
#     U.pop(0)
#     V.pop(0)
#
# plt.show()


# maybe useful to copy:
# for i in range(len(timesteps)):
#     fig,ax=plt.subplots()
#     ax.streamplot(X, Y, U[i], V[i], color=U[i], linewidth=2, cmap='autumn')
#     ax.set_title(f'U V Stream Plot at {timesteps[i]}')
#     plt.show()
#
#     fig, ax = plt.subplots()
#     ax.quiver(X, Y, U[i], V[i])
#     ax.set_title(f'U V Quiver Plot at {timesteps[i]}')
#     plt.show()
