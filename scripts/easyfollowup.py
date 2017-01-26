import numpy as np

# file containing modes found in Step1:
SumFile = '/home/andrea/Desktop/work/DMBayes/chains/step1/1summary.txt'
data1 = np.loadtxt(SumFile)

# file containing the Global Null Evidence found in Step2:
Global = '/home/andrea/Desktop/work/DMBayes/chains/step2/2stats.dat'
data2 = open(Global,'r')

#Null log likelihood : extracted from 2stats.dat
for line in data2.readlines()[:1]:
	if line.startswith('Global Evidence'):
		nullloglike = float(line.strip().split()[2])

dim = 0
for i in range(len(data1)):
	if data1[i,-2] > nullloglike:
		dim += 1

params = np.zeros([dim,7])

n=0
for i in range(len(data1)):
	if data1[i,-2] > nullloglike:
		for j,k in zip(range(7),range(113,120)):
			params[n,j] = data1[i,k]

		n += 1

# Name and path for output follow up file :
followupfilename = '/home/andrea/Desktop/work/DMBayes/Follow_up/Modes_above_NullEvid.txt'
fufile=open(followupfilename,'w')
fufile.write('#File to be read in by Follow-up run \n#Contains modes found by RunI that were identified as PS by model comparison \n#Author: Andrea Chiappo \n \n#Number of PS modes \n')
fufile.write('Number_PS = ' + str(dim) + '\n' + '\n')
fufile.write('#PS params' + '\n')
temp = ' '.join(map(str, params[:,0]))
fufile.write('l = ' + temp + '\n')
temp = ' '.join(map(str, params[:,1]))
fufile.write('b = ' + temp + '\n')
temp = ' '.join(map(str, params[:,2]))
fufile.write('log_N0 = ' + temp + '\n')
temp = ' '.join(map(str, params[:,3]))
fufile.write('E0 = ' + temp + '\n')
temp = ' '.join(map(str, params[:,4]))
fufile.write('alpha_PS = ' + temp + '\n')
temp = ' '.join(map(str, params[:,5]))
fufile.write('beta_PS = ' + temp + '\n')
temp = ' '.join(map(str, params[:,6]))
fufile.write('Inv_Ec = ' + temp + '\n')

fufile.write('#B.F. likelihood value - only important for resume = T, empty otherwise \n#Watch out with the signs. Any actual lnlike value you put in here should be positive \n#BF_Like < 0. means there is no such value \nBF_Like = -1.0')
fufile.close()
