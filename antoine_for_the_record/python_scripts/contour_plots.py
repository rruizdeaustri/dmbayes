#!/usr/bin/python

# This script generates contour plots in (l, b) space
# from .txt files output by calclike.f90 routine contour_plots
# plots for 15 energy bins, for each of model, data, residual

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import math

# Customize path names here
homedir = '/Users/francoisrecanati/Desktop/Stage_Imperial/10PSdat/'
contourpath = homedir + '10PSdat/'
contourpath_true = contourpath + 'pseudotrue/'
savedir = contourpath
outfiles = np.array(['Data.png','Model.png','Residual.png','Likelihood.png'])
prefixes = np.array(['data_','model_','residual_','likelihood_'])


# Number of l,b bins and Energy bins are hard-coded in for now--FIX later
Ebins = 15
bins = 60

# Set arrays with l,b values, containing midpoints of bins
ROI = np.array([[-7.5,-7.5],[7.5,7.5]]) #  [[l_lower,b_lower],[l_upper,b_upper]]
l_int = (ROI[1,0]-ROI[0,0])/bins
b_int = (ROI[1,1]-ROI[0,1])/bins
lvalues = np.arange(ROI[0,0] + b_int/2, ROI[1,0], l_int)
bvalues = np.arange(ROI[0,1] + b_int/2, ROI[1,1], b_int)

# Declare arrays to hold data
plot_pts = np.empty([bins,bins,Ebins,4])
fract_residual = np.empty([bins,bins,Ebins])
likelihood_true = np.empty([bins,bins,Ebins])

# Load in text files and fill data arrays
for j in range(4):
  for i in range(Ebins):
    if i in range(9):
     infilename = contourpath + prefixes[j] + ' ' +  str(i+1) + '.dat'
    else:
     infilename = contourpath + prefixes[j] + str(i+1) + '.dat'
    print(infilename)
    plot_pts[:,:,i,j] = np.loadtxt(infilename)

#Load in likelihood set for true point
for i in range(Ebins):
  if i in range(9):
    infile = contourpath_true + prefixes[3] + ' ' + str(i+1) + '.dat'
  else:
    infile = contourpath_true + prefixes[3] + str(i+1) + '.dat'
  likelihood_true[:,:,i] = np.loadtxt(infile)

data = plot_pts[:,:,:,0]
model = plot_pts[:,:,:,1]
residual = plot_pts[:,:,:,2]
likelihood = plot_pts[:,:,:,3]

# Adding in frational residual counts/bin (i.e. residual/observed)
for i in range(bins):
  for j in range(bins):
    for k in range(Ebins):
      if (plot_pts[i,j,k,0] != 0.0e0):
        fract_residual[i,j,k] = plot_pts[i,j,k,2]/plot_pts[i,j,k,0]
      else:
        fract_residual[i,j,k] = 0.0e0

# Set zero elements in model and data to minimum nonzero value
min_data = np.min(data[np.nonzero(data)])/10.0
min_model = np.min(model[np.nonzero(model)])/10.0
data[np.where(data == 0.0e0)] = min_data
model[np.where(model == 0.0e0)] = min_model
max_data = np.max(data)
max_model = np.max(model)

print 'Data minimum = ' + str(min_data)
print 'Data maximum = ' + str(max_data)
print 'Model minimum = ' + str(min_model)
print 'Model maximum = ' + str(max_model)

# Plot formatting
residlevs = [-800,-600,-400,-280,-240,-200,-160,-120,-80,-40,0,40] #Change, so these are not hardcoded
#Get contour log-levels for data and model plots
U_bound = math.ceil(math.log10(max_data))
L_bound = math.floor(math.log10(min_data))
Nlevs = U_bound + L_bound + 1

# Plotting
for k in range(4):
  fig = plt.figure()

  for j in range(Ebins):

    ax = fig.add_subplot(3,5,j+1)
    cmap = cm.get_cmap('rainbow')
    # Log scale for data and model plots
    if k in (0,1):
      norm = colors.LogNorm(vmin=min_data, vmax=max_data)
      s = ax.contourf(bvalues, lvalues, plot_pts[:,:,j,k], norm=norm, cmap=cmap)
    #elif (k==1):
    #  norm = colors.LogNorm(vmin=min_model, vmax= max_model*10.0)
    #  s = ax.contourf(bvalues, lvalues, plot_pts[:,:,j,k], norm=norm, cmap=cm.get_cmap('rainbow'))
    else:
      s = ax.contourf(bvalues, lvalues, plot_pts[:,:,j,k], cmap=cmap)
  
    ax.set_ylim(lvalues[0], lvalues[bins-1])
    ax.set_xlim(bvalues[0], bvalues[bins-1])
    ax.set_title('Energy bin: ' + str(j+1))
    ax.set_ylabel('l')
    ax.set_xlabel('b')
    for item in ([ax.title,ax.xaxis.label,ax.yaxis.label]):
      item.set_fontsize(6)
    for xitem in (ax.get_xticklabels() + ax.get_yticklabels()):
      xitem.set_fontsize(6)

  #s.set_clim([np.min(plot_pts[:,:,:,k]), np.max(plot_pts[:,:,:,k])])
  cbar_ax=fig.add_axes([0.925,0.1,0.0125,0.8])
  #fig.colorbar(s,cax=cbar_ax)
  if k in (0,1):
    #s.set_clim(min_data, max_data*10.0)
    #ylabs=['1e-1', '1', '1e1', '1e2', '1e3', '1e4', '1e5']
    fig.colorbar(s, cax=cbar_ax, cmap=cm.get_cmap('rainbow'), norm=colors.LogNorm(vmin=min_data,vmax=max_data))
#    if k==0:  
#      ylabs = ['1e-1', '1', '1e1', '1e2', '1e3', '1e4', '1e5']
#    if k==1:
#      ylabs = ['1e-1','1e1','1e3']
    #cbar_ax.set_yticks([1.e-1,1,1e1,1e2,1e3,1e4,1e5])
#    cbar_ax.set_yticklabels(ylabs)
#  if (k==1):
#    s.set_clim(min_model, max_model*10.0)
#    ylabs=['1e-9','1e-7','1e-5','1e-3','1e-1','1e1','1e3','1e5']
#    cbar_ax.set_yticklabels(ylabs)
  else:
    fig.colorbar(s,cax=cbar_ax)
  fig.subplots_adjust(hspace=0.4,wspace=0.4)
  fig.savefig(savedir + outfiles[k], dpi=100, format='png')

# Plotting of fractional residuals
fig = plt.figure()

for j in range(Ebins):
  ax = fig.add_subplot(3,5,j+1)
  s = ax.contourf(bvalues, lvalues, fract_residual[:,:,j])
  
  ax.set_ylim(lvalues[0], lvalues[bins-1])
  ax.set_xlim(bvalues[0], bvalues[bins-1])
  ax.set_title('Energy bin: ' + str(j+1))
  ax.set_ylabel('l')
  ax.set_xlabel('b')
  for item in ([ax.title,ax.xaxis.label,ax.yaxis.label]):
    item.set_fontsize(6)
  for xitem in (ax.get_xticklabels() + ax.get_yticklabels()):
    xitem.set_fontsize(6)

cbar_ax=fig.add_axes([0.925,0.1,0.01,0.8])
fig.colorbar(s,cax=cbar_ax)
fig.subplots_adjust(hspace=0.4,wspace=0.4)
fig.savefig(savedir + 'Fractional_Residual.png', dpi=100, format='png')

#Plotting difference between true and best-fit likelihoods
fig = plt.figure()

for j in range(Ebins):
  ax = fig.add_subplot(3,5,j+1)
  s = ax.contourf(bvalues, lvalues, likelihood_true[:,:,j] - plot_pts[:,:,j,3])
  
  ax.set_ylim(lvalues[0], lvalues[bins-1])
  ax.set_xlim(bvalues[0], bvalues[bins-1])
  ax.set_title('Energy bin: ' + str(j+1))
  ax.set_ylabel('l')
  ax.set_xlabel('b')
  for item in ([ax.title,ax.xaxis.label,ax.yaxis.label]):
    item.set_fontsize(6)
  for xitem in (ax.get_xticklabels() + ax.get_yticklabels()):
    xitem.set_fontsize(6)

cbar_ax=fig.add_axes([0.925,0.1,0.01,0.8])
fig.colorbar(s,cax=cbar_ax)
fig.subplots_adjust(hspace=0.4,wspace=0.4)
fig.savefig(savedir + 'Likelihood_diff.png', dpi=100, format='png')
