import numpy as np
from matplotlib import pyplot as plt


f=open('newdatafile.txt','r')

imax=int(f.readline())
jmax=int(f.readline())
xlength=float(f.readline())
ylength=float(f.readline())
print(imax,jmax,xlength,ylength,sep=" ")

U=f.readline().split('/')[:-1]
U=np.array(U)
U=100*U.astype(np.double).reshape(imax+2,jmax+2)

V=f.readline().split('/')[:-1]
V=np.array(V)
V=100*V.astype(np.double).reshape(imax+2,jmax+2)

P=f.readline().split('/')[:-1]
P=np.array(P)
P=P.astype(np.double).reshape(imax+2,jmax+2)

X=np.arange(0,xlength+1.1*xlength/(imax),xlength/(imax))
Y=np.arange(0,ylength+1.1*ylength/(jmax),ylength/(jmax))
X,Y=np.meshgrid(X,Y)

fig,ax=plt.subplots()
ax.streamplot(X,Y,U,V,color=U, linewidth=2, cmap='autumn')
plt.show()

fig,ax= plt.subplots()
ax.quiver(X,Y,U,V)
ax.set_title('U V Plot')
plt.show()






