#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('eagle_SFH_logger.txt')
rebuild = np.loadtxt('timesteps_32.txt')
time = rebuild[:,3]
redshifts = time[rebuild[:,11]%2==1]

frac = (data[:,7]-data[:,8])/data[:,8]

plotrebuild=True
if plotrebuild:
    SFH_rebuild = data[:,8][rebuild[:,11]%2==1]
    SFH_active_rebuild =  data[:,9][rebuild[:,11]%2==1]
    SFH_inactive_rebuild =  data[:,10][rebuild[:,11]%2==1]
    frac_rebuild = frac[rebuild[:,11]%2==1]

starformed = data[:,4]/np.max(data[:,4])
starformed *= 0.12
    

#plt.xscale('log')
plt.plot(data[:,3],data[:,7],label='SFH inferred')
plt.plot(data[:,3],data[:,8],label='SFH true')
if plotrebuild:
    plt.plot(redshifts,SFH_rebuild,'or',label='rebuild')
plt.ylim(0,160000)
plt.xlim(0,20)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.show()

plt.xscale('log')
plt.plot(data[:,3],frac,label='fractional difference')
if plotrebuild:
    plt.plot(redshifts,frac_rebuild,'or',label='rebuild')
plt.plot(data[:,3],starformed,label='star formed here')
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.show()

plt.xscale('log')
plt.plot(data[:,3],data[:,5],label='SFH active inferred')
plt.plot(data[:,3],data[:,9],label='SFH active true')
if plotrebuild:
    plt.plot(redshifts,SFH_active_rebuild,'or',label='rebuild')
plt.ylim(0,160000)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.show()

plt.xscale('log')
plt.plot(data[:,3],data[:,7]-data[:,5],label='SFH inactive inferred')
plt.plot(data[:,3],data[:,10],label='SFH inactive true')
if plotrebuild:
    plt.plot(redshifts,SFH_inactive_rebuild,'or',label='rebuild')
plt.ylim(0,160000)
plt.legend()
plt.xlabel('Redshift')
plt.ylabel('Star formation density (arbitrary units)')
plt.show()
