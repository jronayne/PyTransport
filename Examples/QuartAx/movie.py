# -*- coding: utf-8 -*-
#################### generate figures ############################################################
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
bet = np.load('data/bet.npy'); alp= np.load('data/alp.npy'); Bztot = np.load('data/alBetBi.npy'); snaps = np.load('data/snaps.npy');
Pz1tot = np.load('data/alBetPz1.npy'); Pz2tot=np.load('data/alBetPz2.npy'); Pz3tot=np.load('data/alBetPz3.npy')
for snap in range(0,np.size(snaps)):
    
    fnlOut = 5.0/6*Bztot[:,:,snap]/(Pz1tot[:,:,snap]*Pz2tot[:,:,snap] + Pz1tot[:,:,snap]*Pz3tot[:,:,snap] + Pz2tot[:,:,snap]*Pz3tot[:,:,snap])
    k =  1543.9416122555126
    kt = 3*k; k1 = kt/2 - bet*kt/2.; k2 = kt/4*(1+alp+bet); k3 = kt/4*(1-alp+bet)    
    DBi = Bztot[:,:,snap]*( k1*k2*k3 )**2
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


# mayavi 3D plot for fnl
    from mayavi import mlab
    mlab.figure(bgcolor=(1.,1.,1.), fgcolor=None, engine=None, size=(600, 600))
    mlab.surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut),vmax=np.nanmax(fnlOut),warp_scale="auto",opacity=1.0, colormap = 'jet')
    #mlab.contour_surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut), warp_scale='auto')
    mlab.view(azimuth=250, elevation=15); mlab.savefig("movie/alpbetFnlMay"+'{0:03}'.format(snap)+".png")
    mlab.close()
    
# mayavi 3D plot for DBi
    from mayavi import mlab
    mlab.figure(bgcolor=(1.,1.,1.), fgcolor=None, engine=None, size=(600, 600))
    mlab.surf(alp,bet,DBi,vmin=np.nanmin(DBi),vmax=np.nanmax(DBi),warp_scale="auto",colormap = 'jet')
    #mlab.contour_surf(alp,bet,fnlOut,vmin=np.nanmin(fnlOut), vmax=np.nanmax(fnlOut), warp_scale='auto')
    mlab.view(azimuth=250, elevation=15); mlab.savefig("movie/alpbetDBiMay"+'{0:03}'.format(snap)+".png")
    mlab.close()
# from command line run:  ffmpeg -framerate 10 -i alpbetFnlMay%03d.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"-c:v libx264 -r 30 -pix_fmt yuv420p out.mp4

    #mlab.show()
