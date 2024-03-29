import numpy as np
import matplotlib.pyplot as plt

f = open('finaldata/liddrivencavity50.txt', 'r')
seconds=100
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
p50 = np.zeros_like(u50)
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
    p50[i] = np.array(f.readline().split('/')[:-1]).astype(np.double).reshape(imax50 + 2, jmax50 + 2)[1:-1, 1:-1]

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
p100 = np.zeros_like(u100)
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
    p100[i] = np.array(g.readline().split('/')[:-1]).astype(np.double).reshape(imax100 + 2, jmax100 + 2)[1:-1, 1:-1]

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
p200 = np.zeros_like(u200)
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
    p200[i] = np.array(h.readline().split('/')[:-1]).astype(np.double).reshape(imax200 + 2, jmax200 + 2)[1:-1, 1:-1]

h.close()

difftopu = u50[1:,:,:] - u100[1:, ::2, ::2]
diffbotu = u100[1:, ::2, ::2] - u200[1:, ::4, ::4]
qvalueu = np.linalg.norm(difftopu, axis=(1, 2)) / np.linalg.norm(diffbotu, axis=(1, 2))

difftopv = v50[1:,:,:] - v100[1:, ::2, ::2]
diffbotv = v100[1:, ::2, ::2] - v200[1:, ::4, ::4]
qvaluev = np.linalg.norm(difftopv, axis=(1, 2)) / np.linalg.norm(diffbotv, axis=(1, 2))
qvalue= np.linalg.norm(np.sqrt(difftopu*difftopu+difftopv*difftopv),axis=(1,2))/np.linalg.norm(np.sqrt(diffbotu*diffbotu+diffbotv*diffbotv),axis=(1,2))

difftopp = p50[1:,:,:] - p100[1:, ::2, ::2]
diffbotp = p100[1:, ::2, ::2] - p200[1:, ::4, ::4]
qvaluep = np.linalg.norm(difftopp, axis=(1, 2)) / np.linalg.norm(diffbotp, axis=(1, 2))


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

fig, ax = plt.subplots(1,figsize=(20,20))
# ax.plot((timesteps_array50 * delt50)[1:], np.log2(qvalueu),color="black",linewidth=5,label='U')
# ax.plot((timesteps_array50 * delt50)[1:], np.log2(qvaluev),color="blue",linewidth=5,label='V')
ax.plot((timesteps_array50*delt50)[1:],np.log2(qvalue),color='#012F5D',linewidth=7,label='U^2+v^2')
ax.plot((timesteps_array50 * delt50)[1:], np.log2(qvaluep),color="#012F5D",linewidth=7,label='P')
ax.hlines(y=1.1449523463098232,xmin=-5,xmax=100,color="black",linestyle=(0, (1, 10)),linewidth=3)
ax.hlines(y=1.0706206606334314,xmin=-5,xmax=100,color="black",linestyle=(0, (1, 10)),linewidth=3)

ax.set_title("Self-Convergence Test")
ax.set_xlabel("Time t in seconds")
ax.set_ylabel("Convergence order p")
ax.set_ylim((1.05,1.3))
ax.set_xlim(-3,103)

plt.show()
plt.savefig('selfconv.pdf',dpi=300)