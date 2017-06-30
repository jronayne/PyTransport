# -*- coding: utf-8 -*-
#################### generate equilateral bispectrum using MTeasyDQuad and MPI ############################################################
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import interpolate 
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
from numpy.random import uniform, seed
from matplotlib.mlab import griddata
#import numpy.ma as ma

# load 3D data
bet = np.load('data/bet.npy'); alp= np.load('data/alp.npy'); Bztot = np.load('data/alBetBi.npy')
Pz1tot = np.load('data/alBetPz1.npy'); Pz2tot=np.load('data/alBetPz2.npy'); Pz3tot=np.load('data/alBetPz3.npy')
fnlOut = 5.0/6*Bztot[:,:,-1]/(Pz1tot[:,:,-1]*Pz2tot[:,:,-1] + Pz1tot[:,:,-1]*Pz3tot[:,:,-1] + Pz2tot[:,:,-1]*Pz3tot[:,:,-1])
k =  2295.8032000107523
kt = 3*k; k1 = kt/2 - bet*kt/2.; k2 = kt/4*(1+alp+bet); k3 = kt/4*(1-alp+bet)    
DBi = Bztot[:,:,-1]*( k1*k2*k3 )**2
#X=np.linspace(-1,1,100); Y=np.linspace(0,1,1F00)
#Y, X = np.meshgrid(X, Y)
#a=np.reshape(alp,(np.size(alp)));b=reshape(bet,(np.size(bet)));f=reshape(fnlOut,np.size(fnlOut));X=reshape(X,np.size(X));Y=reshape(Y,np.size(Y))
#F=griddata(a,b,f, X,Y,interp='linear')


# 3D plot for fnl
#fig1 = plt.figure(1)
#ax = fig1.add_subplot(1, 1, 1, projection='3d')
#surf = ax.plot_surface(alp, bet, fnlOut, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.00, antialiased=False, vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut))    
#ax.grid(True);  ax.set_ylabel(r'$\beta$',fontsize=15);ax.set_xlabel(r'$\alpha$',fontsize=15);  
##ax.set_title('Slice through reduced bispectrum')

#fig1.colorbar(surf, shrink=0.5, aspect=5); plt.savefig("alpBetFnl1.png")

#  projecton plot
fig2 = plt.figure(2,figsize=(10,6),facecolor='white')

ax1 = fig2.add_subplot(1, 1, 1)
cont = ax1.contourf(alp, bet, fnlOut,10, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.0, antialiased=False)#, vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut))    
contLines = ax1.contour(alp, bet, fnlOut, 10, colors='black', linewidth=.5); ax1.clabel(contLines, inline=1, fontsize=10)
ax1.grid(True);   ax1.set_ylabel(r'$\beta$',fontsize=15);ax1.set_xlabel(r'$\alpha$',fontsize=15); 
#ax1.set_title('Slice through reduced bispectrum')
fig2.colorbar(cont, shrink=0.5, aspect=5); plt.savefig("alpBetFnl2.png")

plt.show(fig2); #plt.show(fig1); 

# mayavi 3D plot for fnl
from mayavi import mlab
mlab.figure(bgcolor=(1.,1.,1.), fgcolor=None, engine=None, size=(600, 600))
mlab.surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut),vmax=np.nanmax(fnlOut),warp_scale="auto",opacity=1.0, colormap = 'jet')
#mlab.contour_surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut), warp_scale='auto')
mlab.view(azimuth=250, elevation=15); mlab.savefig("alpbetFnlMay.png")

mlab.show()


# create 3D plot for DBi
#fig3 = plt.figure(3)
#ax = fig3.add_subplot(1, 1, 1, projection='3d')
#surf = ax.plot_surface(alp, bet, DBi, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.00, antialiased=False, vmin=np.nanmin(DBi), vmax=np.nanmax(DBi))    
#ax.grid(True);  ax.set_ylabel(r'$\beta$',fontsize=15);ax.set_xlabel(r'$\alpha$',fontsize=15);  
##ax.set_title('Slice through reduced bispectrum')
#fig3.colorbar(surf, shrink=0.5, aspect=5); plt.savefig("alpBetDBi1.png")

# projection plot for DBi

fig4 = plt.figure(4,figsize=(10,6),facecolor='white')

ax2 = fig4.add_subplot(1, 1, 1)
cont = ax2.contourf(alp, bet, DBi, 10, rstride=1, cstride=1, cmap=cm.jet,linewidth=0.0, antialiased=False, vmin=np.nanmin(DBi), vmax=np.nanmax(DBi) )   
contLines = ax2.contour(alp, bet, DBi, 10, colors='black', linewidth=.5); ax2.clabel(contLines, inline=1, fontsize=10)
ax2.grid(True);   ax2.set_ylabel(r'$\beta$',fontsize=15);ax2.set_xlabel(r'$\alpha$',fontsize=15); 
#ax1.set_title('Slice through reduced bispectrum')
fig4.colorbar(cont, shrink=0.5, aspect=5); plt.savefig("alpBetDBi2.png")

plt.show(fig4); #plt.show(fig3); 

# mayavi 3D plot for fnl
from mayavi import mlab
mlab.figure(bgcolor=(1.,1.,1.), fgcolor=None, engine=None, size=(600, 600))
mlab.surf(alp,bet,DBi,vmin=np.nanmin(DBi),vmax=np.nanmax(DBi),warp_scale="auto",colormap = 'jet')
#mlab.contour_surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut), warp_scale='auto')
mlab.view(azimuth=250, elevation=15); mlab.savefig("alpbetDBiMay.png")

mlab.show()
