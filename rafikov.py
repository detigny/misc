#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 16:12:47 2020

@author: juanito
"""

#The only function imported from FirsUtils is runner, if you load ds directly
#from your snapshot directory modifying the code a little bit, you can delete 
#that line and everything should be fine to be used by anyone :) 
from FirstUtils import *

import yt
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as col
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import *


from SMBH import *

def surface_density(k, resolution=[600,600], w=20, showy=False,gauss=False):
    
    nmbr = "000" + str(k+1)
    if len(str(k+1))==1:
        nmbr = "0000" + str(k+1)
    elif len(str(k+1))==3:
        nmbr = '00'+str(k+1)
    ds = runner(nmbr)

#We will define surf_dens by integrating the density field w/r the z-axis,
#the outcome is a fixed-resolution-buffer that can be readily graphed
    width = (w,'pc')
    surf = ds.proj("density", 'z')
    extent = []
    frb = surf.to_frb(width,resolution)
    
#We may want to do some Gaussian smoothing to the pixels
    if(gauss):
        frb.apply_gauss_beam(nbeam=30, sigma=2.0)
    
    dens = np.flip(np.array(frb['density'].in_units('g/cm**2')),0)
    
    xpos = np.linspace(-w/2,w/2,resolution[0])
    ypos = np.linspace(-w/2,w/2,resolution[1])
    
    if showy:
        ex =  (-w/2, w/2, -w/2,w/2)
        plt.xlabel('$x\:[pc]$')
        plt.ylabel('$y\:[pc]$')
        plt.imshow(dens,extent=ex)
    
        plt.colorbar().set_label(r'Surface Density $\Sigma(x,y)\:\:\: \left[\frac{g}{cm^2}\right]$')
        return 'plotted'
    
    return dens, xpos, ypos

def radial_surf_dens(radius,dens,xpos,ypos,bins=101,showy=False,m=5):
#way of defining and accessing surface density through a radial parametrization
    
    w = np.max(xpos)-np.min(ypos)
    rpc = np.zeros((np.size(xpos),np.size(ypos)))
    for i in range(np.size(xpos)):
        for j in range(np.size(ypos)):
            rpc[i,j] = np.sqrt(xpos[i]**2 + ypos[j]**2) 

    robj = np.zeros(bins-1)
    dens_grid = np.zeros(bins-1)
    R = np.linspace(0,w/2,bins)
    
    for i in range(np.size(R)-1):
        (xi,yi) = np.where((rpc<=R[i+1]) & (rpc>R[i]))
        count = np.size(xi)
        avg_dens = np.sum(dens[xi,yi])/count
        dens_grid[i] = avg_dens
        robj[i] = (R[i+1]+R[i])/2

    step = (bins-1)/4
    step = int(step)
    ministep = int(step/2)+1
        
    r1 = robj[:ministep]
    d1 = dens_grid[:ministep]
    r2 = robj[ministep-1:step]
    d2 = dens_grid[ministep-1:step]
    r3 = robj[step-1:2*step]
    d3 = dens_grid[step-1:2*step]
    r4 = robj[2*step-1:4*step]
    d4 = dens_grid[2*step-1:4*step]
    
    #We do spline and differentiate in every interval
    spl1 = UnivariateSpline(r1, d1, k=m)
    spl2 = UnivariateSpline(r2, d2, k=m)
    spl3 = UnivariateSpline(r3, d3, k=m)
    spl4 = UnivariateSpline(r4, d4, k=m)
    sg = np.concatenate((spl1(r1[:-1]),spl2(r2[:-1]),spl3(r3[:-1]),spl4(r4)))

    if showy:
        plt.figure()
        plt.xlabel(r'$R\:\:\:[pc]$')
        plt.ylabel(r'$\Sigma(R)\:\:\: \left[\frac{g}{cm^2}\right]$')
        
        plt.plot(robj,dens_grid,'.')
        # plt.plot(r1,spl1(r1))
        # plt.plot(r2,spl2(r2))
        # plt.plot(r3,spl3(r3))
        # plt.plot(r4,spl4(r4))
        plt.plot(robj,sg)
        plt.show()
        return 'plotted'

    if (0 < R) | (R < r1[-1]):
        return spl1(R)
    elif (r2[0] <= R) | (R < r2[-1]):
        return spl2(R)
    elif (r3[0] <= R) | (R < r3[-1]):
        return spl3(R)
    elif (r4[0] <= R) | (R <= r4[-1]):
        return spl4(R)
    else:
        return 'Seems like the radius is outside our range defined at the FRB'
        
    

def omega(k, resolution=[600,600], w=20, showy=False,gauss=False):
    
    nmbr = "000" + str(k+1)
    if len(str(k+1))==1:
        nmbr = "0000" + str(k+1)
    elif len(str(k+1))==3:
        nmbr = '00'+str(k+1)
    ds = runner(nmbr)
    
#Same rationale as with density, but here we integrate tangential velocity 
#with a weighting by the density field and then normalizing. Omega is then
#calculated in norm, by division of the radius
    width = (w,'pc')
    surf = ds.proj("tangential_velocity", 'z',weight_field='density')
    extent = []
    frb = surf.to_frb(width,resolution)
    if(gauss):
        frb.apply_gauss_beam(nbeam=30, sigma=2.0)
    
    vels = np.flip(np.array(frb['tangential_velocity'].in_units('cm/s')),0)
    
    xpos = np.linspace(-w/2,w/2,resolution[0])
    ypos = np.linspace(-w/2,w/2,resolution[1])
    rpc = np.zeros((resolution[0],resolution[1]))
    for i in range(np.size(xpos)):
        for j in range(np.size(ypos)):
            rpc[i,j] = np.sqrt(xpos[i]**2 + ypos[j]**2)
    rcm = rpc*3.086e18
    
    if showy:
        ex =  (-w/2, w/2, -w/2,w/2)
        plt.xlabel('$x\:[pc]$')
        plt.ylabel('$y\:[pc]$')
        plt.imshow(vels/rcm,extent=ex,norm=col.LogNorm())
        plt.colorbar().set_label(r'$\Omega(R)\:\:\: \left[\frac{1}{s}\right]$')
        return 'plotted'

    return vels/rcm, xpos, ypos, rpc

def setup_omega(k, bins=101, resolution=[600,600], w=20, m=5,showy=False):
    oms, x, y, r = omega(k,resolution,w)
    
    robj = np.zeros(bins-1)
    omegs = np.zeros(bins-1)
    
    #Rotational speed radial profile
    R = np.linspace(0,w/2,bins)
    for i in range(np.size(R)-1):
        (xi,yi) = np.where((r<=R[i+1]) & (r>R[i]))
        count = np.size(xi)
        avg_rot = np.sum(oms[(xi,yi)])/count
        omegs[i] = avg_rot
        robj[i] = (R[i+1]+R[i])/2
    
    step = (bins-1)/4
    step = int(step)
    ministep = int(step/2)+1

    
    #4-interval-split for interpolating different regimes
    
    r1 = robj[1:ministep]
    om1 = omegs[1:ministep]
    r2 = robj[ministep-1:step]
    om2 = omegs[ministep-1:step]
    r3 = robj[step-1:2*step]
    om3 = omegs[step-1:2*step]
    r4 = robj[2*step-1:4*step]
    om4 = omegs[2*step-1:4*step]
    
    #We do a spline interpolation and differentiate in every interval
    spl1 = UnivariateSpline(r1, om1, k=m)
    spl_der1 = spl1.derivative()
    spl2 = UnivariateSpline(r2, om2, k=m)
    spl_der2 = spl2.derivative()
    spl3 = UnivariateSpline(r3, om3, k=m)
    spl_der3 = spl3.derivative()    
    spl4 = UnivariateSpline(r4, om4, k=m)
    spl_der4 = spl4.derivative()
    
    #We generate derivative data points and we interpolate in intervals again
    dp1 = spl_der1(r1)
    dp2 = spl_der2(r2)
    dp3 = spl_der3(r3)
    dp4 = spl_der4(r4)
    dps = np.concatenate((dp1[:-1],dp2[:-1],dp3[:-1],dp4))
    
    if showy:
        plt.figure()
        plt.xlabel(r'$R\:\:\:[pc]$')
        plt.ylabel(r'$\Omega(R)\:\:\: \left[\frac{1}{s}\right]$')
        
        plt.plot(robj,omegs,'.')
        plt.plot(r1,spl1(r1))
        plt.plot(r2,spl2(r2))
        plt.plot(r3,spl3(r3))
        plt.plot(r4,spl4(r4))
        plt.show()
    
    deriv1 = UnivariateSpline(r1,dp1)
    deriv2 = UnivariateSpline(r2,dp2)
    deriv3 = UnivariateSpline(r3,dp3)
    deriv4 = UnivariateSpline(r4,dp4)
    
    deriv = np.concatenate((deriv1(r1[:-1]),deriv2(r2[:-1]),deriv3(r3[:-1]),deriv4(r4)))
    
    if showy:    
        plt.figure()
        plt.xlabel(r'$R\:\:\:[pc]$')
        plt.ylabel(r'$\frac{d\Omega(R)}{dR}\:\:\: \left[\frac{1}{s}\right]$')
        plt.plot(robj[1:],dps,'.')
        plt.plot(robj[1:],deriv)
        plt.show()
        return 'Plotted some data :)'
    
    return deriv1, deriv2, deriv3, deriv4, r1,r2,r3,r4

def r_der_omega(k, R, bins=101, resolution=[600,600], w=20, preload=False):
#If you want to evaluate in several R's, preload the derivative from setup
#Else it's VERY slow to load the Omega-frb from the setup everytime  
    
    if preload==False:
        d1,d2,d3,d4,r1,r2,r3,r4 = setup_omega(k,bins,resolution,w)
    else:
        d1,d2,d3,d4,r1,r2,r3,r4 = preload
    
    if (0 < R) | (R < r1[-1]):
        return d1(R)
    elif (r2[0] <= R) | (R < r2[-1]):
        return d2(R)
    elif (r3[0] <= R) | (R < r3[-1]):
        return d3(R)
    elif (r4[0] <= R) | (R <= r4[-1]):
        return d4(R)
    else:
        return 'Seems like the radius is outside our range defined at the FRB'
    
def momentum_flux(k,R,ders,b):
#The whole point here is to profile this, so it's always with preload :)
    pass


