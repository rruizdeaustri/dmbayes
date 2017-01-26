# indexing an array of points allpts, assuming
# that its columns (vector indices) are:
# vals, civ, gen, ..., ptsrc

# column numbers after calling load_pts are
# 0: values, 1:civilization, 2:generation, 3: l, 4: b, 
# 5: log(N0), 6:E0, 7:alpha, 8:beta, 9:Inv Ec

import numpy as np
import matplotlib.pyplot as plt
import math

def true_pts(filename):
  pts = np.loadtxt(filename, usecols=(0,1,2,3,4,5,6)) 
  return pts

# read the important points out of the .sam file
def load_pts(filename):
    #filename = 'BG_3D_PS_7D_NP' + str(NP) + '.sam'
    pts = np.loadtxt(filename, usecols=(1,2,3,35,36,37,38,39,40,41,42))
    return pts
        
# returns points for a given civ, generation range, ptsrc, and parameter list
def usepts(allpts, civ, genlo, genhi, ptsrc, vec_index): 
    usepts = allpts[np.where(allpts[:,1] == civ)]
    usepts = usepts[np.where(usepts[:,2] <= genhi)]
    usepts = usepts[np.where(usepts[:,2] >= genlo)]
    usepts = usepts[np.where(usepts[:,-1] == ptsrc)]
    usepts = usepts[np.where(usepts[:,0] <= 1.0e+31)] # added this line to remove infinite points for plotting over all generations
    usepts = usepts[:, vec_index]
    return usepts

# Plot the evolution of the effective population size for each point source with generation
def plot_NPeff(allpts, Nps, civ, genlo, genhi, savedir='plots/'):
	
	plt.figure()
	plt.title('Civilization ' + str(civ))
	plt.xlabel('Generation')
	plt.ylabel('Effective Population Size') 	

	pts = np.empty([Nps+1,genhi])	
	civpts = allpts[np.where(allpts[:,1] == civ)]

	for gen in range(genlo,genhi+1):
		genpts = civpts[np.where(civpts[:,2] == gen)] 
		pts[0,gen-1]=gen
		for ptsrc in range(Nps):
    			#PSpts = genpts[np.where(genpts[:,10] == ptsrc)]
			pts[ptsrc+1,gen-1]=np.size(np.where(genpts[:,10] == ptsrc+1))

	for ps in range(Nps):
		plt.plot(pts[0,:], pts[ps+1,:], label='PS ' + str(ps+1))

	plt.legend(loc=2)
	plt.savefig(savedir + 'civ' + str(civ) + '_NP.png', dpi=100, format='png')

#Find the best fit parameter vector for each point source 
def find_min(allpts,civ,genlo,genhi,ptsrc):

        pts=usepts(allpts,civ,genlo,genhi,ptsrc,(0,2,3,4,5,6,7,8,9))

        minval = min(pts[:,0])
        minlocs = np.where(pts[:,0] == minval)

        gen = np.take(pts[:,1], minlocs[0])
        lmin = np.take(pts[:,2], minlocs[0])
        bmin = np.take(pts[:,3], minlocs[0])
        Nmin = np.take(pts[:,4], minlocs[0])
        Emin = np.take(pts[:,5], minlocs[0])
        alphamin = np.take(pts[:,6], minlocs[0])
        betamin = np.take(pts[:,7], minlocs[0])
        Ecmin = np.take(pts[:,8], minlocs[0])
        print('-------------------------------------------------------')
	print('Point source ' + str(ptsrc) + ':')
#        print( 'Minimum value = ' + str(minval) + ' at (' +
#                str(lmin[0]) + ', ' + str(bmin[0])  + ', ' + str(Nmin[0])  + ', ' + str(Emin[0])
#                + ', ' + str(alphamin[0]) + ', ' + str(betamin[0])  + ', ' + str(Ecmin[0])
#                + ')' )
#        print( 'for Generation ' + str(gen[0]))
        print( 'Best-fit Likelihood = ' + str(minval))
        print( 'at {' + str(Nmin[0]) + ', ' + str(Emin[0]) + ', ' + str(alphamin[0]) + ', ' + str(betamin[0]) + ', ' + str(Ecmin[0]) + '}')
        print( '(l, b) = ' + '(' + str(lmin[0]) + ', '  + str(bmin[0]) + ')')

# Return best-fit parameter array {logN0, E0, alpha, beta, 1/Ec}
def bestfit_pt(allpts,civ,genlo,genhi,ptsrc):

	pts = usepts(allpts,civ,genlo,genhi,ptsrc,(0,2,3,4,5,6,7,8,9))

	minval = min(pts[:,0]) 
	minlocs = np.where(pts[:,0] == minval)

        lmin = np.take(pts[:,2], minlocs[0])
        bmin = np.take(pts[:,3], minlocs[0])
	N0min = np.take(pts[:,4], minlocs[0])
	E0min = np.take(pts[:,5], minlocs[0])
	alphamin = np.take(pts[:,6], minlocs[0])
	betamin = np.take(pts[:,7], minlocs[0])
	inv_Ecmin = np.take(pts[:,8], minlocs[0])

	bestfit_pt=np.array([lmin[0], bmin[0], N0min[0], E0min[0], alphamin[0], betamin[0], inv_Ecmin[0], minval])

	return bestfit_pt

# Plot the evolution of the likelihood value 
def plot_evol_like(allpts, Nps, civ, genlo, genhi, savedir='plots/'):
 
        # Added minimum best-fit values
        minvals = np.empty([Nps])
        for ptsrc in range(Nps):	
         pts=usepts(allpts,civ,genlo,genhi,ptsrc+1,(0,2,3,4,5,6,7,8,9))
         minvals[ptsrc] = min(pts[:,0])

        generations = np.arange(genlo, genhi+1)
        npts = np.size(generations)
        likepts = np.empty([npts,Nps])

        # Added values of BF and true likelihood, needs to be done manually
        true_like = 113006.602476925*np.ones(npts)
        BF_like = 114888.581743668*np.ones(npts)

	civpts = allpts[np.where(allpts[:,1] == civ)]
	#civpts = civpts[np.where(civpts[:,0] <= 1.0e+31)]

        for i in range(npts):
          genpts = civpts[np.where(civpts[:,2] == generations[i])]
          for ptsrc in range(Nps):
#            print('Generation ' + str(generations[i]) + ', Point source ' + str(ptsrc+1))
            PSpts = genpts[np.where(genpts[:,10] == ptsrc+1)]
            likepts[i, ptsrc] = min(PSpts[:,0])

        plt.figure()
        plt.xlabel('Generation')
        plt.ylabel(r'$\Delta$ BF log like')
	for ps in range(Nps):
          plt.plot(generations, likepts[:, ps]-minvals[ps], label='PS ' +str(ps+1))
      
        plt.yscale('log')
	plt.xlim(genlo,genhi)
	plt.ylim(0.01,1.e6)
	plt.legend(loc=1)
	plt.savefig(savedir + 'BFLike_diff.png', dpi=100, format='png')

        plt.figure()
        plt.xlabel('Generation')
        plt.ylabel('Best-fit Likelihood')
        for ps in range(Nps):
          plt.plot(generations, likepts[:,ps], label='PS ' + str(ps+1))
        plt.plot(generations, true_like, 'm--', label='True Point')
        #plt.plot(generations, BF_like, 'c--', label='Best-fit Point')
        plt.ylim(1.e5, 5.e5)
        plt.legend(loc=1)
        plt.savefig(savedir + 'BFLike.png', dpi=100, format='png')

# Plot the fractional improvement of the best-fit likelihood value for a given point source
def plot_like_improve(allpts,ptsrc,civ,genlo,genhi,savedir='plots/'):

	plt.figure()
	plt.xlabel('Generation')

        npts = genhi-genlo+1

	likepts = np.empty([2,npts])
	likeimprove = np.empty([2,npts-1])
	civpts = allpts[np.where(allpts[:,1] == civ)]
	#civpts = civpts[np.where(civpts[:,0] <= 1.0e+10)]

	for gen in range(npts):
	  genpts = civpts[np.where(civpts[:,2] == gen + genlo)]
	  likepts[0,gen] = gen + genlo
	  PSpts = genpts[np.where(genpts[:,10] == ptsrc)]
          likepts[1,gen] = min(PSpts[:,0])
        
        for gen in range(npts-1):
          likeimprove[0,gen] = gen + genlo + 1 
          likeimprove[1,gen] = (likepts[1,gen] - likepts[1,gen+1])/likepts[1,gen]
        
        PSpts = likeimprove[:,likeimprove[1,:] != 0.0]
             
        plt.plot(PSpts[0,:], PSpts[1,:], 'b.', label='PS' + str(1))

        plt.yscale('log')
	plt.xlim(genlo+1, genhi)
	plt.legend(loc=1)
	plt.savefig(savedir+ 'Likelihood_Improvement.png', dpi=100, format='png') 
        

# Plot the evolution of the mean likelihood value
# Also superimposed best-fit likelihoods for each point source 
def mean_like_evol(allpts, Nps, civ, genlo, genhi, savedir='plots/'):
	
	plt.figure()
	plt.xlabel('Generation')
	plt.ylabel('Fractional Likelihood Improvement') 	

        generations = np.arange(genlo, genhi+1)
        npts = np.size(generations)
	meanlike = np.empty([npts])
        PSlike = np.empty([npts,Nps])
        fracmean = np.empty([npts-1])
        fracPS = np.empty([npts-1,Nps])	 

	civpts = allpts[np.where(allpts[:,1] == civ)]
	#civpts = civpts[np.where(civpts[:,0] <= 1.0e+10)]

	for i in range(npts):
	  genpts = civpts[np.where(civpts[:,2] == generations[i])] 
	  meanlike[i] = np.mean(genpts[:,0])
          for ptsrc in range(Nps):
	    PSpts = genpts[np.where(genpts[:,10] == ptsrc+1)]
            PSlike[i,ptsrc] = min(PSpts[:,0])

        for j in range(npts-1):
          fracmean[j] = (meanlike[j] - meanlike[j+1])/meanlike[j]
          for ptsrc in range(Nps):
            fracPS[j,ptsrc] = (PSlike[j,ptsrc] - PSlike[j+1,ptsrc])/PSlike[j,ptsrc]

	plt.plot(generations[1:npts], fracmean, 'm.', label='Mean Likelihood')
        for ps in range(Nps):
          plt.plot(generations[1:npts], fracPS[:,ps], linecolors[ps], label= 'PS ' +str(ps+1))

        plt.yscale('log')
	plt.xlim(genlo,genhi)
	plt.legend(loc=1)
        plt.savefig(savedir + 'LikeEvol.png', dpi=100, format='png')

#Plot evolution of mean likelihood value
def plot_mean_like(allpts, civ, genlo, genhi, savedir='plots/'):
	
	plt.figure()
	plt.xlabel('Generation')
	plt.ylabel('Mean Likelihood') 	

        generations = np.arange(genlo, genhi+1)
        npts = np.size(generations)
	meanlike = np.empty([npts])

	civpts = allpts[np.where(allpts[:,1] == civ)]
	#civpts = civpts[np.where(civpts[:,0] <= 1.0e+10)]

	for i in range(npts):
	  genpts = civpts[np.where(civpts[:,2] == generations[i])] 
	  meanlike[i] = np.mean(genpts[:,0])

	plt.plot(generations, meanlike, 'm', label='Mean Likelihood')

        plt.yscale('log')
	plt.xlim(genlo,genhi)
	#plt.legend(loc=1)
        plt.savefig(savedir + 'MeanLike.png', dpi=100, format='png')

# Plot of proposed convergence criterion
def plot_cnvg(allpts,Nps,civ,genlo,genhi,savedir='plots/'):

	plt.figure()
	plt.xlabel('Generation')
        plt.ylabel('Fractional Mean Likelihood Improvement')

        linecolours = np.array(['k.','m.']) 
        pltlabels = np.array(['weighted', 'unweighted'])

        w = 10
        generations = np.arange(genlo, genhi+1)
        npts = np.size(generations)
        meanpts = np.empty([npts]) 
	BFlikepts = np.empty([npts,Nps])
        fracpts = np.empty([npts-1,Nps+1])
        Nbins = int(np.floor(npts/w))
        binpts = np.empty([Nbins])
	#plotpts = np.empty([Nbins,5])
        plotpts = np.empty([npts-10,Nps+2])

        for k in range(Nbins):
          binpts[k] = generations[(k+1)*w]
          
	civpts = allpts[np.where(allpts[:,1] == civ)]
	#civpts = civpts[np.where(civpts[:,0] <= 1.0e+10)]

        for i in range(npts):
          genpts = civpts[np.where(civpts[:,2] == generations[i])]
          meanpts[i] = np.mean(genpts[:,0])
          for ps in range(Nps):
            PSpts = genpts[np.where(genpts[:,10] == ps+1)]
            BFlikepts[i,ps] = min(PSpts[:,0])

        for m in range(npts-1):
          for ps in range(Nps):
            fracpts[m,ps] = (BFlikepts[m,ps] - BFlikepts[m+1,ps])/BFlikepts[m,ps]
          fracpts[m,Nps] = (meanpts[m] - meanpts[m+1])/meanpts[m]

        for j in range(npts-10):
          jL = j 
          jR = j + w - 1
          for ps in range(Nps):
            plotpts[j,ps] = generations[jR+1]*np.mean(fracpts[jL:jR+1,ps])
          plotpts[j,Nps] = generations[jR+1]*np.mean(fracpts[jL:jR+1,Nps])
          plotpts[j,Nps+1] = np.mean(fracpts[jL:jR+1,Nps])
        #for j in range(Nbins):
        #  gL = j*w
        #  gR = (j+1)*w
        #  print binpts[j]
        #  for ps in range(3):
        #    plotpts[j,ps] = binpts[j]*np.mean(fracpts[gL:gR,ps]) 
        #  plotpts[j,3] = binpts[j]*np.mean(fracpts[gL:gR,3])
        #  plotpts[j,4] = np.mean(fracpts[gL:gR,3])
 
        for col in (Nps,Nps+1):
          plt.plot(generations[10:], plotpts[:,col], linecolours[col], label=pltlabels[col])

        plt.yscale('log')
	plt.xlim(genlo+1, genhi)
	plt.legend(loc=1)
	plt.savefig(savedir+ 'Convergence_Criteria.png', dpi=100, format='png') 

# scatterplot of points in given civ, gen range 
# Modified from param_scplot() to give single figure with plots over all parameter combinations and point sources
def scplot(allpts, truevals, ptsrc, civ, genlo, genhi, save=False, savedir='plots/'):

   # determine the names of parameter & label plot
    paramnames = {3:'l', 4:'b', 5:'log(N0)', 6:'E0', 
                       7:'alpha', 8:'beta', 9:'Inv_Ec'}
 
    paramindices = np.array([3,5,7,8])

    fig=plt.figure()

    #for ptsrc in (1,2,3):
    for i in (1,2,3,4): 

        param1=paramindices[i-1]
	param2=paramindices[i-1]+1
        thesepts = usepts(allpts, civ, genlo, genhi, ptsrc, (0, param1, param2))

        ax=fig.add_subplot(2,2,i)
    	s=ax.scatter(thesepts[:,1], thesepts[:,2], c=thesepts[:,0], s=8, edgecolors='none', alpha=0.8)
    #Plot true point
        ax.scatter(truevals[ptsrc-1,param1-3], truevals[ptsrc-1,param2-3], s=120, c='y', marker='^', edgecolors='r')
    # find best point, print to screen, and plot
    # minlocs is an array (best point occurs in successive generations)
        minval = min(thesepts[:,0])
        minlocs = np.where(thesepts[:,0] == minval)
        xmin = np.take(thesepts[:,1], minlocs[0])
        ymin = np.take(thesepts[:,2], minlocs[0])
        ax.scatter(xmin[0], ymin[0], s=120, c='r', marker='^',edgecolors='cyan')

        if param1 in paramnames:
          param1name=paramnames[param1]
          ax.set_xlabel(param1name)
        else:
          param1name = ''
        
        if param2 in paramnames:
          param2name=paramnames[param2]
          ax.set_ylabel(param2name)
        else:
          param2name = ''
        
        ax.set_title('PS ' + str(ptsrc) + ' (' + param1name + ', ' + param2name + ')')

        #for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        #  item.set_fontsize(8)
	
	#for xitem in (ax.get_xticklabels()+ ax.get_yticklabels()):
        #  xitem.set_fontsize(6)
       

    fig.subplots_adjust(hspace=0.4,wspace=0.4)
    cbar_ax=fig.add_axes([0.925,0.1,0.025,0.8])
    for citem in (cbar_ax.get_yticklabels() + cbar_ax.get_xticklabels()):
      citem.set_fontsize(8)
    fig.colorbar(s,cax=cbar_ax)
	
    # save figure
    if save:
        fig.savefig(savedir + 'PS' + str(ptsrc) + '_params'  + '.png', dpi=100, format='png')

#Plot the evolution of the best fit parameters with generation
def param_evol(allpts, civ, genlo, genhi, ptsrc, savedir='plots/'):

   # determine the names of parameter & label plot
    paramnames = {3:'l', 4:'b', 5:'log(N0)', 6:'E0', 
                       7:'alpha', 8:'beta', 9:'Inv_Ec', 0:'Likelihood'}
    paramindices = np.array([3,4,5,6,7,8,9,0])

    linecolours = np.array([['b','b--'],['g','g--'],['r','r--']])

    fig=plt.figure()
    
    generations = np.arange(genlo, genhi+1)
    npts = np.size(generations)
    BF_params = np.empty([npts,8])

    #Constant array containing true points
    trueparams = np.ones([npts,7])
    for i in range(7):
      trueparams[:,i] = truevals[ptsrc-1,i]*trueparams[:,i]

    for i in range(npts):
      BF_params[i,:] = bestfit_pt(allpts, civ, genlo, generations[i], ptsrc)

    for j in range(7):
      ax=fig.add_subplot(4,2,j+1)
      ax.plot(generations[:], BF_params[:,j], linecolours[ptsrc-1,0], trueparams[:,j], linecolours[ptsrc-1,1])
      paramlabel=paramnames[paramindices[j]]
      ax.set_ylabel(paramlabel)
      ax.set_xlabel('Generation')
  
      #ax.set_title()

      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
        item.set_fontsize(8)	
      for xitem in (ax.get_xticklabels()+ ax.get_yticklabels()):
        xitem.set_fontsize(6)
    
    ax=fig.add_subplot(4,2,8)
    ax.plot(generations[:], BF_params[:,7], linecolours[ptsrc-1,0])
    ax.set_ylabel(paramnames[paramindices[7]])
    ax.set_xlabel('Generations')   
    ax.set_yscale('log')  
 
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label]):
      item.set_fontsize(8)	
    for xitem in (ax.get_xticklabels()+ ax.get_yticklabels()):
      xitem.set_fontsize(6)
    
    fig.subplots_adjust(hspace=0.4,wspace=0.4)
    fig.savefig(savedir + 'BFparam_evol_PS' + str(ptsrc)  + '.png', dpi=100, format='png')

# Plot energy spectra of best-fit point and true point
def plot_spectra(allpts, truevals, ptsrc, civ, genlo, genhi, save=False, savedir='/plots/'):

	E = np.arange(1000,100000,10)
	BF_params=np.empty([8])
	dNdE_BF = np.empty([np.size(E)])
	dNdE_true = np.empty([np.size(E)])
	
	BF_params = bestfit_pt(allpts,civ,genlo,genhi,ptsrc)
          #print('PS' + str(ptsrc))
          #print(BF_params[ptsrc-1,:])
        for i in range(np.size(E)):
	  dNdE_BF[ i] = dNdE(E[i], BF_params[2], BF_params[3], BF_params[4], BF_params[5], BF_params[6])
	  dNdE_true[i] = dNdE(E[i], truevals[2], truevals[,3], truevals[4], truevals[5], truevals[6])
	
        plt.figure()
	plt.xlabel('Energy (MeV)')

	linecolors = np.array([['b','g','r','m','k'],['b--','g--','r--','m--','k--']])

        plt.plot(E, dNdE_BF[:], 'b', label='PS ' + str(ptsrc))
        plt.plot(E, dNdE_true[:], 'b--')

	plt.xscale('log')
        plt.yscale('log')
	plt.xlim(1000,10000)
	plt.ylim(1.e-13,1.e-8)
        plt.legend(loc=1)
	plt.savefig(savedir  + 'Energy_Spectra_PS' + str(ptsrc) + '.png', dpi=100, format='png')


# Function to calculate energy spectrum from point source spectral parameters
# Note that N0 here is log(N0) in the functional definition
def dNdE(E, N0, E0, alpha, beta, inv_Ec):

	dNdE = 10.0**(N0)*((E/E0)**(-alpha - beta*math.log10(E/E0)))*math.exp(-(E-E0)*inv_Ec)	
        return dNdE
