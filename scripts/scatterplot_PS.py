import numpy as np
import math as m
import matplotlib.pyplot as plt
import sys

SumFile = '/home/andrea/Desktop/work/DMBayes/chains/step1/1summary.txt'
TxtFile = '/home/andrea/Desktop/work/plot/46PScatalogtrueparam.txt'

#path to .txt file :
data=np.loadtxt(SumFile)
trupts=np.loadtxt(TxtFile, usecols=(1,2,5))
fig=plt.figure()
ax=fig.add_subplot(1,1,1)
plt.xlabel('l [deg]')
plt.ylabel('b [deg]')
s = ax.scatter(data[:,113],data[:,114],s=50,c=data[:,-1],label='MultiNest')
t = ax.scatter(trupts[:,0],trupts[:,1],s=50, c = 'grey', marker = 'D', edgecolors = 'pink',label='True PS')
cbar_ax=fig.add_axes([0.9,0.1,0.025,0.8])
cb = fig.colorbar(s,cax=cbar_ax)
cb.set_label(r'$\chi^2$',rotation=0)
ax.legend(scatterpoints=1,loc='upper right',fontsize=12)
fig.suptitle('6D Multinest [nCdims = 2 , live = 3000 , tol = 0.01 , eff = 0.3]')
fig.savefig('/home/andrea/Desktop/work/plot/MultiNest/MNmodes_46TruePS' + '.png', dpi=100, format='png')
plt.show()
