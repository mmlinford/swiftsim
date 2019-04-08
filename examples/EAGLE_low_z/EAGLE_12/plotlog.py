#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import sys

case = int(sys.argv[1])
cosmo=True
if case ==0:
    # Nocosmo
    data = np.loadtxt('eagle_nocosmo_SFH_logger.txt')
    rebuild = np.loadtxt('timesteps_18.txt')
    name = 'nocosmo'
    cosmo=False
elif case ==1:
    # Regular
    data = np.loadtxt('eagle_SFH_logger.txt')
    rebuild = np.loadtxt('timesteps_16.txt')
    name = 'regular'
else:
    # Force rebuild
    data = np.loadtxt('eagle_forcerebuild_SFH_logger.txt')
    rebuild = np.loadtxt('timesteps_22.txt')
    name = 'forcerebuild'

redshiftarray = data[:,3]
if cosmo==False:
    redshiftarray = data[:,1]
redshifts = redshiftarray[rebuild[:,11]%2==1]
plotrebuild=True

frac = (data[:,7]-data[:,8])/data[:,8]
if plotrebuild:
    SFH_rebuild = data[:,8][rebuild[:,11]%2==1]
    SFH_active_rebuild =  data[:,9][rebuild[:,11]%2==1]
    SFH_inactive_rebuild =  data[:,10][rebuild[:,11]%2==1]
    frac_rebuild = frac[rebuild[:,11]%2==1]
    

if cosmo:
    plt.xscale('log')
plt.plot(redshiftarray,data[:,7],label='SFH inferred')
plt.plot(redshiftarray,data[:,8],label='SFH true')
if plotrebuild:
    plt.plot(redshifts,SFH_rebuild,'or',label='rebuild')
plt.ylim(0,4000)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.savefig('./plots/Total_SFH_'+name+'.png')
plt.show()

converted = data[:,4]/np.max(data[:,4])

if cosmo:
    plt.xscale('log')
plt.plot(redshiftarray,frac,label='SFH inferred')
plt.plot(redshiftarray,converted*0.0014,label='star formed')
if plotrebuild:
    plt.plot(redshifts,frac_rebuild,'or',label='rebuild')
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.savefig('./plots/deviation_SFH_'+name+'.png')
plt.show()

if cosmo:
    plt.xscale('log')
plt.plot(redshiftarray,data[:,5],label='SFH active inferred')
plt.plot(redshiftarray,data[:,9],label='SFH active true')
plt.plot(redshiftarray,converted*np.max(data[:,9]),label='star formed')
if plotrebuild:
    plt.plot(redshifts,SFH_active_rebuild,'or',label='rebuild')
plt.ylim(0,4000)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.savefig('./plots/Active_SFH_'+name+'.png')
plt.show()

if cosmo:
    plt.xscale('log')
plt.plot(redshiftarray,data[:,7]-data[:,5],label='SFH inactive inferred')
plt.plot(redshiftarray,data[:,10],label='SFH inactive true')
if plotrebuild:
    plt.plot(redshifts,SFH_inactive_rebuild,'or',label='rebuild')
plt.ylim(0,4000)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.savefig('./plots/Inactive_SFH_'+name+'.png')
plt.show()
