import numpy as np
import matplotlib.pyplot as plt

f = open('finaldata/liddrivencavity50.txt', 'r')
seconds=50
skipsteps=100
delt=0.001

imax50 = int(f.readline())
jmax50 = int(f.readline())
xlength50 = float(f.readline())
ylength50 = float(f.readline())
delt50 = float(f.readline())
print(imax50,jmax50,xlength50,ylength50,delt50,sep=" ")
delt50=delt

timesteps50 = int(seconds / delt50 / skipsteps) + 1
print(timesteps50)
u50 = np.zeros((timesteps50, imax50, jmax50), dtype=np.double)
v50 = np.zeros_like(u50)
timesteps_array50 = np.zeros((timesteps50))
for i in range(4):
    f.readline()

for i in range(timesteps50):
    line = f.readline()
    if not line:
        break
    timesteps_array50[i] = int(line)
    u50[i] = np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax50 + 2, jmax50 + 2)[1:-1, 1:-1]
    v50[i] = np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax50 + 2, jmax50 + 2)[1:-1, 1:-1]
    f.readline()

f.close()

g = open('finaldata/liddrivencavity100.txt', 'r')

imax100 = int(g.readline())
jmax100 = int(g.readline())
xlength100 = float(g.readline())
ylength100 = float(g.readline())
delt100 = float(g.readline())
delt100=delt
timesteps100 = int(seconds / delt100/ skipsteps) + 1
print(timesteps100)
u100 = np.zeros((timesteps100, imax100, jmax100), dtype=np.double)
v100 = np.zeros_like(u100)
timesteps_array100 = np.zeros((timesteps100))
for i in range(4):
    g.readline()

for i in range(timesteps100):
    line = g.readline()
    if not line:
        break
    timesteps_array100[i] = int(line)
    u100[i] = np.array(g.readline().split('/')[:-1]).astype(np.double).reshape(imax100 + 2, jmax100 + 2)[1:-1, 1:-1]
    v100[i] = np.array(g.readline().split('/')[:-1]).astype(np.double).reshape(imax100 + 2, jmax100 + 2)[1:-1, 1:-1]
    g.readline()

g.close()

h = open('finaldata/liddrivencavity200.txt', 'r')

imax200 = int(h.readline())
jmax200 = int(h.readline())
xlength200 = float(h.readline())
ylength200 = float(h.readline())
delt200 = float(h.readline())
delt200=delt
timesteps200 = int(seconds / delt200/ skipsteps) + 1
print(timesteps200)
u200 = np.zeros((timesteps200, imax200, jmax200), dtype=np.double)
v200 = np.zeros_like(u200)
timesteps_array200 = np.zeros((timesteps200))
for i in range(4):
    h.readline()

for i in range(timesteps200):
    line = h.readline()
    if not line:
        break
    timesteps_array200[i] = int(line)
    u200[i] = np.array(h.readline().split('/')[:-1]).astype(np.double).reshape(imax200 + 2, jmax200 + 2)[1:-1, 1:-1]
    v200[i] = np.array(h.readline().split('/')[:-1]).astype(np.double).reshape(imax200 + 2, jmax200 + 2)[1:-1, 1:-1]
    h.readline()

h.close()

difftop = u50[1:,:,:] - u100[1:, ::2, ::2]
diffbot = u100[1:, ::2, ::2] - u200[1:, ::4, ::4]
qvalue = np.linalg.norm(difftop, axis=(1, 2)) / np.linalg.norm(diffbot, axis=(1, 2))

fig, ax = plt.subplots()
ax.plot((timesteps_array50 * delt50)[1:], qvalue)
plt.show()