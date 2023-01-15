import numpy as np
from matplotlib import pyplot as plt

f = open('datafiles/newdatafile.txt', 'r')

# read parameters

imax = int(f.readline())
jmax = int(f.readline())
xlength = float(f.readline())
ylength = float(f.readline())
print(imax, jmax, xlength, ylength, sep=" ")

# initial fields

U = f.readline().split('/')[:-1]
U = np.array(U)
U = 100 * U.astype(np.double).reshape(imax + 2, jmax + 2)
U = np.swapaxes(U, 0, 1)

V = f.readline().split('/')[:-1]
V = np.array(V)
V = 100 * V.astype(np.double).reshape(imax + 2, jmax + 2)
V = np.swapaxes(V, 0, 1)

P = f.readline().split('/')[:-1]
P = np.array(P)
P = P.astype(np.double).reshape(imax + 2, jmax + 2)
P = np.swapaxes(P, 0, 1)

# geometry

FLAG = f.readline().split('/')[:-1]
FLAG = np.array(FLAG).astype(np.double)
malen = [1. if x == 16 else 0. for x in FLAG]
malen = np.array(malen).reshape(imax + 2, jmax + 2).astype(np.double)
malen = np.swapaxes(malen, 0, 1)
konvertieren = np.zeros((imax + 2, jmax + 2, 3))
konvertieren[malen > 0.5] = [1, 1, 1]
konvertieren[malen < 0.5] = [0, 0, 0]
konvertieren=konvertieren[::-1,:]

malfig, malax = plt.subplots()
malax.imshow(konvertieren, interpolation='nearest')
malax.set_title("Geometry")
plt.show()


X = np.arange(0, xlength + 1.1 * xlength / imax, xlength / imax)
Y = np.arange(0, ylength + 1.1 * ylength / jmax, ylength / jmax)
X, Y = np.meshgrid(X, Y)

# fig, ax = plt.subplots()
# ax.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
# plt.show()
#
# fig, ax = plt.subplots()
# ax.quiver(X, Y, U, V)
# ax.set_title('U V Plot')
# plt.show()

# now read additional data
ind = 2  # jeder 4. datenpunkte wird geplottet
plot_number = 4
while True:

    line = f.readline()
    if not line:
        break
    timestep = int(line)
    U = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    V = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    P = np.swapaxes(np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax + 2, jmax + 2), 0, 1)
    if plot_number % 1 == 0:
        fig, ax = plt.subplots()
        ax.streamplot(X, Y, U, V, color=U, linewidth=2, cmap='autumn')
        ax.set_title(f'U V Stream Plot at {timestep}')
        plt.show()

        fig, ax = plt.subplots()
        ax.quiver(X[::ind, ::ind], Y[::ind, ::ind], U[::ind, ::ind], V[::ind, ::ind])
        ax.set_title(f'U V Quiver Plot at {timestep}')
        plt.show()
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