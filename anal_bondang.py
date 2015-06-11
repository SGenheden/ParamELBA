# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to plot bond and angle distributions from a trajectory

The trajectory is assumed to be created by md_bondang.py

The force field definition, i.e. the bond and angles to be tracked
are stored in an XML format and read by the classes in the elbalib
module.

A reference structure need to be given and this is used to find CG
molecules defined in the XML file. 

Examples
--------
anal_bondang.py -f bondang_elba.txt bondang_slipid.txt  -l elba slipid -r popc_elba.pdb -x ~/Dropbox/Code/ParamELBA/elba.xml
"""

import argparse
import sys

import matplotlib.pylab as plt
import MDAnalysis as md
import numpy as np
import scipy.stats as stats
import scipy.optimize as opt

import elbalib
import colors

def _fit_dists(distlist,labels,fflist) :
    """
    Fit the distributions and compute 
    """
    
    def _harmonic(xdata, k, eq):
        return k*(xdata-eq)*(xdata-eq)
    
    nseries = distlist[0].shape[1]
    RT      = 303.0*8.314/4.184/1000.0

    for iseries in range(nseries):
        
        fftype = fflist[iseries].type
        isbond = isinstance(fftype,elbalib.BondType)
        if isbond:
            print "ff: k=%.3f r_eq=%.3f"%(fftype.k*fftype.kf,fftype.r0),
        else:
            print "ff: k=%.3f theta_eq=%.3f"%(fftype.k*fftype.kf,fftype.theta0),

        for i,(dists,label) in enumerate(zip(distlist,labels)):
            l       = float(dists.shape[0])
            if isbond:
                y = dists[:,iseries]
            else:
                y = np.cos(np.deg2rad(dists[:,iseries]))
            h,e = np.histogram(y,bins=l/40.0)
            kde = stats.gaussian_kde(y)
            kde.covariance_factor = lambda : 0.25
            kde._compute_covariance()
            reffunc = -RT*np.log(kde(e))
            p0 = [RT/y.var(),y.mean()]
            popt,pcov = opt.curve_fit(_harmonic,e,reffunc,p0)
            if isbond:
                print " %s: k=%.3f r_eq=%.3f"%(label,popt[0],popt[1]),
            else:
                print " %s: k=%.3f theta_eq=%.3f"%(label,popt[0],np.rad2deg(np.arccos(popt[1]))),
        print ""

def _plot_dists(distlist,labels,fflist,figure,trans=None):
    """
    Plot the distribution of a series in subplots
    """
    nseries = distlist[0].shape[1]
    ncols   = 3.0
    nrows   = np.ceil(nseries / ncols)
    RT      = 303.0*8.314/4.184/1000.0

    for iseries in range(nseries):
        a = figure.add_subplot(nrows,ncols,iseries+1)
        minval = 1000
        maxval = -1000
        stdval = 0
        for i,(dists,label) in enumerate(zip(distlist,labels)):
            l       = float(dists.shape[0])
            if trans is None:
                y = dists[:,iseries]
            else:
                y = trans(dists[:,iseries])
            h,e = np.histogram(y,bins=l/40.0)
            density = stats.gaussian_kde(y)
            density.covariance_factor = lambda : 0.25
            density._compute_covariance()
            print y.std(),
            a.plot(e,density(e),color=colors.color(i),label=label)
            stdval = max(stdval,np.round(2.0*y.std(),0),0.5)
            minval = min(minval,np.round(np.floor(e.min()),0))
            maxval = max(maxval,np.round(np.ceil(e.max()),0))
        mean,std = fflist[iseries].type.statmoments(RT)
        print mean,std
        minval = min(minval,np.floor(mean-2.0*std))
        maxval = max(maxval,np.ceil(mean+2.0*std))
        stdval = max(stdval,np.round(2.0*std,0),0.5)
        x = np.arange(minval,maxval+stdval,stdval)
        a.plot(x,fflist[iseries].type.distribution(x,RT),color=colors.color(len(labels)),label="ff")
        if iseries == 0 : a.legend(loc=1,fontsize=8)
        a.set_yticks([])
        a.set_xticks(x)
    figure.tight_layout()

def _read_series(filename,nbonds):
    """
    Read series of data output by md_bondang.py
    """
    data = []
    with open(filename,"r") as f:
        line = f.readline()
        while line:   
            data.append(line.strip().split())
            line = f.readline()
    data = np.array(data,dtype=float)
    bonddata = data[:,1:nbonds+1]
    angdata  = data[:,nbonds+1:]
    return bonddata,angdata

if __name__ == '__main__':

    print " ".join(sys.argv)

    # Command-line input
    parser = argparse.ArgumentParser(description="Analyse distribution of bonds and angles")
    parser.add_argument('-f', '--file', nargs="+", help="the trajectory of the bond and angles")
    parser.add_argument('-l', '--labels', nargs="+", help="the labels of the different trajectories")
    parser.add_argument('-r', '--ref', help="a CG reference file", default="ref.pdb")
    parser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    parser.add_argument('-o', '--out', help="the output prefx", default="bondang")
    args = parser.parse_args()

    # Load the force field definitions
    ff = elbalib.Elba()
    ff.load(args.xml)

    # Load the reference CG universe
    refuni = md.Universe(args.ref)

    # Find which residue the bond and angles were tracked
    ffmol = None
    for r in refuni.selectAtoms("all").residues:
        ffmol = ff.find_molecule(r.name)
        if ffmol is not None and ffmol.bonds and ffmol.angles : break
    if ffmol is None : raise Exception("Could not find any of the FF molecules in the reference structure")
    nbonds = len(ffmol.bonds)

    # Now we can read data, assume it was created by md_bondang.py
    anglist = []
    bondlist = []
    for filename in args.file:
        bd,ad = _read_series(filename,nbonds)
        bondlist.append(bd)
        anglist.append(ad)

    # Create one plot for bonds and one for angles
    f = plt.figure(1)
    _plot_dists(bondlist,args.labels,ffmol.bonds,f)
    f.savefig(args.out+"_bonddist.png",format="png")
    print ""    

    f = plt.figure(2)
    def _trans_ang(y):
        return np.cos(np.deg2rad(y))
    _plot_dists(anglist,args.labels,ffmol.angles,f,trans=_trans_ang)
    f.savefig(args.out+"_angdist.png",format="png")  
    print ""    

    # Estimate ff parameters
    _fit_dists(bondlist,args.labels,ffmol.bonds)  
    print ""
    _fit_dists(anglist,args.labels,ffmol.angles)  
