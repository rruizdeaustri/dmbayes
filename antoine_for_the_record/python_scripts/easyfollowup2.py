#!/usr/bin/python
import re
import numpy as np
#.stats file to read in the modes :
statsfile=open('/Users/francoisrecanati/Desktop/PS3D_BG3Dstats.dat','r')
# Name and path for output follow up file :
followupfilename = '/Users/francoisrecanati/Desktop/Follow_ups/test2.txt'
def read_in_num(string):
	for line in statsfile:
		if string in line:
			number=int(str(re.findall(r'[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?',line)[0]))
			return number

num_modes= read_in_num('Total Modes Found')
loglike=np.empty(num_modes)

#Null log likelihood : to be written in by hand !
nullloglike=-0.7792282e4

ii=0
for line in statsfile:
	if 'Strictly Local Evidence' in line:
		loglike[ii]= float(line[26:51])
		ii=ii+1
#print(loglike)


accept=[]
for i in range(num_modes):
	if loglike[i] > nullloglike:
		accept.append(int(i))
num_PS=len(accept)

k_vals=np.empty(num_PS)
params=np.empty([num_PS,7])
#Lines=['']*num_modes #np.empty(num_modes)
#Lines=np.empty(num_modes)
for i in  range(num_PS):
	k_vals[i] = 2 + 3*accept[i]
k=0
#print('k vals = ' + str(k_vals))
statsfile.close()
statsfile=open('/Users/francoisrecanati/Desktop/PS3D_BG3Dstats.dat','r')
for line in statsfile:
	if '  32  ' in line:
		k=k+1
#		print(str(line))
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,0]=float(line[7:32])
	if '  33  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,1]=float(line[-24:])
	if '  34  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,2]=float(line[-24:])
	if '  35  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,3]=float(line[-24:])
	if '  36  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,4]=float(line[-24:])
	if '  37  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,5]=float(line[-24:])
	if '  38  ' in line:
		for i in range(num_PS):
			if (k == k_vals[i]):
				params[i,6]=float(line[-24:])
				
print('k = ' + str(k))
#print(params)	

print(params[:,0])
statsfile.close()

fufile=open(followupfilename,'w')

#fufile=open(followupfilename,'w')
fufile.write('#File to be read in by Follow-up run \n#Contains modes found by RunI that were identified as PS by model comparison \n#Author: Charlotte Strege \n \n#Number of PS modes \n')
fufile.write('Number_PS = ' + str(num_PS) + '\n' + '\n')
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


