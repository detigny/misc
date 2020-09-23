#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 10:56:26 2020

@author: juanito
"""

import yt
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate as intr
from numpy.linalg import norm as norm
from scipy.interpolate import griddata


class SMBHpair(object):
    
    def __init__(self, color='orange', simulation_n=1, dt=1, orbital=150,comp_path='none',seps='none',title='none'):
        """Object initialization of course!"""
        
        self.accretions = 'none'
        self.positions = 'none'
        self.binarymass = 'none'
        self.eddington_rate = 'none'
        self.binary_instance_info = 'No snapshot query yet'
        self.merger_alert = np.nan

        if type(simulation_n)==int:
            self.N = simulation_n
        elif type(simulation_n)==list:
            self.N = simulation_n[-1]-simulation_n[0]
            
        self.title = ''
        if title != 'none':
            self.title = title+': '
        
        #For convenience we allow setting separations manually
        self.separations = seps
        
        #dt is the time-step size in the simulation run in 10[kyr] units
        #orbital is the orbital time of the binary in [kyr] units
        self.torb = orbital
        self.dt = dt
        if type(simulation_n)==int:
            self.times = np.linspace(0,(simulation_n-1)*dt*10/orbital,simulation_n)
        elif type(simulation_n)==list:
            self.times = np.linspace((simulation_n[0]-1)*dt*10/orbital,
                                     (simulation_n[-1]-1)*dt*10/orbital,(simulation_n[-1]-simulation_n[0]))

        if comp_path == 'none':
        #This assumes I'm getting info from my HD in a specified format :D
            root = '//Volumes/juanHD1/data_'+color+'/output_'
        elif comp_path != 'none':
        #If I want to use another directory for simulations - specify full path
            root = comp_path
        
        self.sink_directories = []
        if type(simulation_n)==int:
            for i in range(simulation_n):
                number = '000' + str(i+1)
                if len(str(i+1))==1:
                    number = "0000" + str(i+1)
                elif len(str(i+1))==3:
                    number = "00" + str(i+1)
                self.sink_directories.append(root+number+'/sink_'+number+'.csv')
         
        elif type(simulation_n)==list:
            for i in range(simulation_n[0],simulation_n[-1]):
                number = '000' + str(i+1)
                if len(str(i+1))==1:
                    number = "0000" + str(i+1)
                elif len(str(i+1))==3:
                    number = "00" + str(i+1)
                self.sink_directories.append(root+number+'/sink_'+number+'.csv')            

            
    def spatial(self, a_graph=False, p_graph=False, path_graph=False):
        #Here we set de binary separations, with a toggle option for positions
        
        simulation_n = self.N
        seps = np.zeros(simulation_n)
        pos = np.zeros((simulation_n,2,3))
        for i in range(simulation_n):
            sinkinfo = np.genfromtxt(self.sink_directories[i],delimiter=',')
            
            try:
                np.shape(sinkinfo)[1]
                
                pos[i][0][:] = np.array([sinkinfo[0,2],sinkinfo[0,3],sinkinfo[0,4]])
                pos[i][1][:] = np.array([sinkinfo[1,2],sinkinfo[1,3],sinkinfo[1,4]])
                seps[i] = norm((pos[i][0]-pos[i][1]))
                    
            except:
                if np.isnan(self.merger_alert):
                    self.merger_alert = i+1
                    print("There was a binary merger at snapshot " + str(self.merger_alert) + "!")
                    print("Ask for more info avalaible with the .merged routine")
                    seps[i]=0
                    pos[i][0][:] = np.array([sinkinfo[2],sinkinfo[3],sinkinfo[4]])
                    pos[i][1][:] = np.nan
                                
        if self.separations=='none':
            self.separations = seps
        seps = self.separations
        self.positions = pos
        
        
        if a_graph == True:                
            x = self.times
            plt.figure()
            plt.xscale('log')
            plt.yscale('log')
            plt.ylim([0.01,1.5])
            plt.xlim([0.1,100])
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel('$a_{bin}/a_0$')
            plt.title(self.title+'Binary separation')
            plt.plot(x,np.array(seps)/seps[0],linewidth=0.8)
            plt.show()

        if p_graph == True:                
            x = self.times
            y1 = pos[:,0,0]
            y2 = pos[:,1,0]
            plt.figure()
            plt.xscale('log')
            plt.xlim([0.1,10])
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel('x [pc]')
            plt.title(self.title+'X coordinate evolution')
            plt.plot(x,y1,linewidth=0.8)
            plt.plot(x,y2,linewidth=0.8)
            plt.show()
            
            x = self.times
            z1 = pos[:,0,2]
            z2 = pos[:,1,2]
            plt.figure()
            plt.xscale('log')
            plt.xlim([0.1,10])
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel('z [pc]')
            plt.title(self.title+'Z coordinate evolution')
            plt.plot(x,z1,linewidth=0.8)
            plt.plot(x,z2,linewidth=0.8)
            plt.show()
            
            y1 = pos[:,0,1]
            y2 = pos[:,1,1]
            plt.figure()
            plt.xscale('log')
            plt.xlim([0.1,10])
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel('y [pc]')
            plt.title(self.title+'Y coordinate evolution')
            plt.plot(x,y1,linewidth=0.8)
            plt.plot(x,y2,linewidth=0.8)
            plt.show()
            
        if path_graph == True:                
            x = pos[:][0][0]
            y = pos[:][0][1]
            plt.figure()
            plt.xscale('log')
            plt.xlabel('x [pc]')
            plt.ylabel('y [pc]')
            plt.title(self.title+'X-Y path SMBH-1')
            plt.plot(x,y,linewidth=0.8)
            plt.show()

            x = pos[:][1][0]
            y = pos[:][1][1]
            plt.figure()
            plt.xscale('log')
            plt.xlabel('x [pc]')
            plt.ylabel('y [pc]')
            plt.title(self.title+'X-Y path SMBH-1')
            plt.plot(x,y,linewidth=0.8)
            plt.show()
            
        pass
    
    
    
    def feeding(self, m_graph=False, acc_graph=False, scaled_accretion=False):
        #Here we see the individual SMBH masses, accretion and Edd accretion
        
        if scaled_accretion:
            print("You asked for accretion to be given in Eddington scale")
        
        sinkmass = np.zeros((self.N,2))
        sinkaccretions = np.zeros((self.N,2))
        eddington = np.zeros((self.N,2))
        
        condition = True
        if self.N==1:
            print("No accretion because there is no time step to compare!")
            condition = False
        
        sinkaccretions[0][0] = np.nan
        sinkaccretions[0][1] = np.nan
        
        for i in range(self.N):
            sinkinfo = np.genfromtxt(self.sink_directories[i],delimiter=',')
            try:
                np.shape(sinkinfo)[1]
                
                sinkmass[i][0] = sinkinfo[0,-2]
                sinkmass[i][1] = sinkinfo[1,-2]
                if condition and i>0:
                    sinkaccretions[i][0] = (sinkmass[i][0]- sinkmass[i-1][0])/(self.dt*10*1000)
                    sinkaccretions[i][1] = (sinkmass[i][1]- sinkmass[i-1][1])/(self.dt*10*1000)
                
            except:
                if np.isnan(self.merger_alert):
                    self.merger_alert = i+1
                    print("There was a binary merged at snapshot" + str(self.merger_alert) + "!")
                    print("Ask for more info avalaible with the .merger routine")
                    
                    sinkmass[i][0] = sinkinfo[-2]
                    sinkmass[i][1] = np.nan
                    eddington[i][0] = sinkinfo[-2]/45e6
                    eddington[i][1] = np.nan
                if condition and i>=self.merger_alert:
                    sinkaccretions[i][0] = (sinkmass[i][0]- sinkmass[i-1][0])/(self.dt*10*1000)
                    sinkaccretions[i][0] = np.nan
            
            if np.isnan(self.merger_alert)==False:
                i = self.merger_alert - 1
                sinkaccretions[i][0] = (sinkmass[i][0]- sinkmass[i-1][0] - sinkmass[i-1][1])/(self.dt*10*1000)
            
            if scaled_accretion:
                sinkaccretions[i][0] = sinkaccretions[i][0]/(10*sinkmass[i][0]/45e6)
                sinkaccretions[i][1] = sinkaccretions[i][1]/(10*sinkmass[i][1]/45e6)
                
            self.accretions = sinkaccretions
            self.binarymass = sinkmass
            
        if m_graph == True:                
            y1 = np.array(self.binarymass[:,0])/self.binarymass[0,0]
            y2 = np.array(self.binarymass[:,1])/self.binarymass[0,1]
            x = self.times
            plt.figure()
   #         plt.xscale('log')
    #        plt.xlim([0.1,10])
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel('$M_{BH}/M_0$')
            plt.title(self.title+'Mass evolution of the binary')
            plt.plot(x,y1,linewidth=0.8,color='blue')
            plt.plot(x,y2,linewidth=0.8,color='red')
            plt.show()
                
        if acc_graph == True:
            y1 = np.array(self.accretions[:,0])
            y2 = np.array(self.accretions[:,1])
            x = self.times
 #           plt.xlim([0.1,10])
            titulo = "Average accretion in $10^6M_\odot/yr$ units"
            ylab = "$\dot{M}\:[10^6M_\odot/yr]$"
            if scaled_accretion:
                titulo = titulo = "Average accretion in Eddington units"
                ylab = "$\dot{M}/\dot{M}_{Edd}$"
                
            plt.figure()
  #          plt.xscale('log')
            plt.xlabel('$t/t_{orb}$')
            plt.ylabel(ylab)
            plt.title(titulo)
            plt.plot(x,y1,linewidth=0.8,color='blue')
            plt.plot(x,y2,linewidth=0.8,color='red')
            plt.show()
                
            pass

           
            
    def velocities(self):
        return 1
             
    def merged(self):
        #Merger info and stats!
        if np.isnan(self.merger_alert):
            return "Sadly, I don't have any information to give you about mergers"
        mtime = str((self.merger_alert-1)*self.dt*10)
        print("The merger happened at " + mtime + " kyrs after simulation init")
        
        try:
            posit = str(self.positions[self.merger_alert-1][0])
            print("The position where the merger happened is: " + posit)
            sepa = str(self.separations[self.merger_alert-2][0])
            print("The separation of the SMBHs before merging was: " + sepa)
        except:
            "No merger position available. Maybe try calling .spatial"
        
        try:
            mmass = str(self.binarymass[self.merger_alert-1][0])
            print("The mass of the merged SMBHs turns to be " + mmass)
        except:
            "No merger mass available. Maybe try calling .feeding"
        pass
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        