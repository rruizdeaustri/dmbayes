#!/usr/bin/python

# Python sript for PS plots, using vis_pts.py

import ptsrc_plt_Nps as PS

civ=1
genlo=2
genhi=300
plotdir='plots/10PS/Try1/'

pts=PS.load_pts('/homes/tapajo/dupuisg/DMBayes/Work/PS_test/de/10PS/BG_3D_PS_7D_try1.sam')
truevals=PS.true_pts('/homes/tapajo/dupuisg/DMBayes/True_10PS.txt')

for ptsrc in range(10):
  PS.scplot(pts, truevals, ptsrc+1, civ, genlo,genhi, save=True, savedir=plotdir)

PS.plot_spectra(pts, truevals, civ, genlo, genhi, save=True, savedir=plotdir)

