#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

Data = np.loadtxt('output_SFH_logger.txt')
rebuild = np.loadtxt('timesteps_16.txt')
redshifts=rebuild[:,3][rebuild[:,11]%2==1]

plt.plot(Data[:,0],Data[:,7]*1.022690e-02,label='Total SFH')
plt.plot(Data[:,0],Data[:,8]*1.022690e-02,label='Total SFH')
plt.plot(Data[:,0],Data[:,5]*1.022690e-02,label='Active SFH')
plt.plot(Data[:,0],Data[:,7]*1.022690e-02-Data[:,5]*1.022690e-02,label='inactive SFH')
for i in range(0,len(redshifts)):
    plt.axvline(x=redshifts[i],color='r')
plt.xlabel('Time (z)')
plt.ylabel('SFH [$\\rm M_\odot \\rm yr^{-1}$]')
#plt.ylim(0,50)
plt.legend()
plt.savefig('SFH_lowz_eagle12.png')
plt.show()
plt.close()


plt.plot(Data[:,0],Data[:,7]*1.022690e-02,label='Total SFH logger')
plt.plot(Data[:,0],Data[:,8]*1.022690e-02,label='Total SFH true')
plt.xlabel('Time (z)')
plt.ylabel('SFH [$\\rm M_\odot \\rm yr^{-1}$]')
plt.legend()
plt.show()

plt.plot(Data[:,0],Data[:,5]*1.022690e-02,label='Active SFH logger')
plt.plot(Data[:,0],Data[:,9]*1.022690e-02,label='Active SFH true')
plt.xlabel('Time (z)')
plt.ylabel('SFH [$\\rm M_\odot \\rm yr^{-1}$]')
plt.legend()
plt.show()

plt.plot(Data[:,0],Data[:,7]*1.022690e-02-Data[:,5]*1.022690e-02,label='Inactive SFH logger')
plt.plot(Data[:,0],Data[:,10]*1.022690e-02,label='Inactive SFH true')
plt.xlabel('Time (z)')
plt.ylabel('SFH [$\\rm M_\odot \\rm yr^{-1}$]')
plt.legend()
plt.show()

CSFH_Mstar = np.cumsum(Data[:,4])
CSFH_SFRdt = np.cumsum(Data[:,6])

plt.plot(Data[:,0],CSFH_Mstar*1e10,label='CSFH stars converted')
plt.plot(Data[:,0],CSFH_SFRdt*1e10,label='CSFH SFR tracers')
plt.xlabel('Time (z)')
plt.ylabel('Total stellar mass [$\\rm M_\odot$]')
plt.legend()
plt.savefig('CSFH_lowz_eagle12.png')
plt.show()
plt.close()
