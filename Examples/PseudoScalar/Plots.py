#################### generate equilateral bispectrum using MTeasyDQuad and MPI ############################################################
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import numpy.ma as ma
# make up some randomly distributed data

bet = np.load('data/bet.npy'); alp= np.load('data/alp.npy'); Bztot = np.load('data/alBetBi.npy')
Pz1tot = np.load('data/alBetPz1.npy'); Pz2tot=np.load('data/alBetPz2.npy'); Pz3tot=np.load('data/alBetPz3.npy')
fnlOut = 5.0/6*Bztot[:,:,-1]/(Pz1tot[:,:,-1]*Pz2tot[:,:,-1] + Pz1tot[:,:,-1]*Pz3tot[:,:,-1] + Pz2tot[:,:,-1]*Pz3tot[:,:,-1])

fig1 = plt.figure(1)

ax = fig1.add_subplot(1, 1, 1, projection='3d')
surf = ax.plot_surface(alp, bet, fnlOut, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.00, antialiased=False, vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut))    
ax.grid(True);  ax.set_ylabel(r'$\beta$',fontsize=15);ax.set_xlabel(r'$\alpha$',fontsize=15);  
#ax.set_title('Slice through reduced bispectrum')

fig1.colorbar(surf, shrink=0.5, aspect=5); plt.savefig("alpBetEq2.png")


fig2 = plt.figure(2,figsize=(10,6),facecolor='white')

ax1 = fig2.add_subplot(1, 1, 1)
cont = ax1.contourf(alp, bet, fnlOut, 10, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.0, antialiased=False, vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut))    
contLines = ax1.contour(alp, bet, fnlOut, 10, colors='black', linewidth=.5); ax1.clabel(contLines, inline=1, fontsize=10)
ax1.grid(True);   ax1.set_ylabel(r'$\beta$',fontsize=15);ax1.set_xlabel(r'$\alpha$',fontsize=15); 
#ax1.set_title('Slice through reduced bispectrum')
fig2.colorbar(cont, shrink=0.5, aspect=5); plt.savefig("alpBetEq3.png")

plt.show(fig1); plt.show(fig2); 

from mayavi import mlab
mlab.figure(bgcolor=(1.,1.,1.), fgcolor=None, engine=None, size=(600, 600))
mlab.surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut),vmax=np.nanmax(fnlOut),warp_scale="auto",colormap = 'jet')
#mlab.contour_surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut), warp_scale='auto')
mlab.view(azimuth=250, elevation=15)
mlab.show()

mlab.savefig("test.png")

data = loadtxt('dataT/EqBi.dat')

fig4 = plt.figure(4)
kOut = data[0,:]; fnlOut = data[1,:]; Bztot = data[2,:]; Pztot = data[3,:]
plt.plot(np.log(kOut/kOut[0]), fnlOut, linewidth=2)
#title('Reduced bispectrum in equilateral configuration',fontsize=15) 
grid(True); plt.legend(fontsize=15); ylabel('$fNL$', fontsize=20)
xlabel(r'$\log(k/k_{\rm pivot})$', fontsize=15); grid(True); plt.legend(fontsize=15); plt.xlim(min(np.log(kOut/kOut[0])),max(np.log(kOut/kOut[0]))); plt.savefig("BiEq.png")
plt.show(fig4)
