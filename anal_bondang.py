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
anal_bondang.py -f bondang_elba.txt bondang_slipid.txt  -l elba slipid -r popc_elba.pdb -x elba.xml
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


class TimeSeries(object):
    """
    Simple class to hold a time series and processes thereof
    """

    def __init__(self, data, transf=None):
        self.data0 = data
        if transf is None:
            self.data = data
        else:
            self.data = transf(data)
        self.len = self.data.shape[0]
        self.mean = self.data.mean()
        self.std = self.data.std()
        self.var = self.data.var()
        self.kde = stats.gaussian_kde(self.data)
        self.kde.covariance_factor = lambda: 0.25
        self.kde._compute_covariance()
        self.histogram, self.edges = np.histogram(self.data, bins=self.len / 40.0)

    def density(self, points=None):
        """
        Return the density of the Kernal density estimator
        """
        if points is None:
            return self.kde(self.edges)
        else:
            return self.kde(points)


def _fit_ff(serieslist, labels, fflist):
    """
    Fit the distributions and compute 
    """

    # This is the function we will be fitting to
    def _harmonic(xdata, k, eq, C):
        return k * (xdata - eq) * (xdata - eq) + C

    nseries = len(fflist)
    RT = 303.0 * 8.314 / 4.184 / 1000.0

    for iseries in range(nseries):
        # Print out the force field parameters
        fftype = fflist[iseries].type
        isbond = isinstance(fftype, elbalib.BondType)
        if isbond:
            print "%s | ff: k=%9.3f r_eq=%9.3f" % (fflist[iseries].atomstring(), fftype.k * fftype.kf, fftype.r0),
        else:
            print "%s | ff: k=%9.3f theta_eq=%9.3f" % (fflist[iseries].atomstring(), fftype.k * fftype.kf, fftype.theta0),
        # Now loop over each of the simulated systems and fit to an harmonic system
        for i, (series, label) in enumerate(zip(serieslist, labels)):
            s = series[iseries]
            reffunc = -RT * np.log(s.density())
            p0 = [RT / s.var, s.mean, s.std * np.sqrt(np.pi * 2.0)]
            popt, pcov = opt.curve_fit(_harmonic, s.edges, reffunc, p0)
            if isbond:
                print " %s: k=%9.3f r_eq=%9.3f" % (label, popt[0], popt[1]),
            else:
                print " %s: k=%9.3f theta_eq=%9.3f" % (label, popt[0], np.rad2deg(np.arccos(popt[1]))),
        print ""


def _plot_dists(serieslist, labels, fflist, figure):
    """
    Plot the distribution of a series in subplots
    """
    nseries = len(fflist)
    ncols = 3.0
    nrows = np.ceil(nseries / ncols)
    RT = 303.0 * 8.314 / 4.184 / 1000.0

    for iseries in range(nseries):
        a = figure.add_subplot(nrows, ncols, iseries + 1)
        minval = 1000
        maxval = -1000
        stdval = 0
        # Plot the density for each simulated system
        for i, (series, label) in enumerate(zip(serieslist, labels)):
            s = series[iseries]
            a.plot(s.edges, s.density(), color=colors.color(i), label=label)         
            stdval = max(stdval, np.round(2.0 * s.std, 0), 0.5)
            minval = min(minval, np.floor(s.edges.min()))
            maxval = max(maxval, np.ceil(s.edges.max()))
        # Calculate the range of the x axis
        mean, std = fflist[iseries].type.statmoments(RT)
        minval = min(minval, np.floor(mean - 2.0 * std))
        maxval = max(maxval, np.ceil(mean + 2.0 * std))
        if isinstance(fflist[iseries].type, elbalib.AngleType):
            minval = max(-1.0,minval)
            maxval = min(1.0,maxval)
        stdval = max(stdval, np.round(2.0 * std, 0), 0.5)
        # Plot the force field ideal distribution
        x = np.arange(minval, maxval + stdval, stdval)
        a.plot(x, fflist[iseries].type.distribution(x, RT), color=colors.color(len(labels)), label="ff")  
        # Add legend and ticks
        if iseries == 0:
            a.legend(loc=1, fontsize=8)
        a.set_yticks([])
        a.set_xticks(x)
        if isinstance(fflist[iseries].type, elbalib.AngleType):
            a.set_xticklabels(np.rad2deg(np.arccos(x)))
            a.invert_xaxis()
    figure.tight_layout()


def _read_series(filename, nbonds):
    """
    Read series of data output by md_bondang.py
    """
    data = []
    with open(filename, "r") as f:
        line = f.readline()
        while line:
            data.append(line.strip().split())
            line = f.readline()
    data = np.array(data, dtype=float)
    bonddata = data[:, 1:nbonds + 1]
    angdata = data[:, nbonds + 1:]
    return bonddata, angdata


def _trans_ang(y,inv=False):
    if not inv:
        return np.cos(np.deg2rad(y))
    else:
        return np.rad2deg(np.arccos(y))

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
        if ffmol is not None and ffmol.bonds and ffmol.angles:
            break
    if ffmol is None:
        raise Exception("Could not find any of the FF molecules in the reference structure")
    nbonds = len(ffmol.bonds)
    nangles = len(ffmol.angles)

    # Now we can read data, assume it was created by md_bondang.py
    anglist = []
    bondlist = []
    for filename in args.file:
        bd, ad = _read_series(filename, nbonds)
        # This transforms that numpy arrays into TimeSeries objects and thereby estimate densities etc
        bondlist.append([TimeSeries(bd[:, i]) for i in range(nbonds)])
        anglist.append([TimeSeries(ad[:, i], _trans_ang) for i in range(nangles)])

    # Create one plot for bonds and one for angles
    f = plt.figure(1)
    _plot_dists(bondlist, args.labels, ffmol.bonds, f)
    f.savefig(args.out + "_bonddist.png", format="png")

    f = plt.figure(2)
    _plot_dists(anglist, args.labels, ffmol.angles, f)
    f.savefig(args.out + "_angdist.png", format="png")

    # Estimate ff parameters
    _fit_ff(bondlist, args.labels, ffmol.bonds)
    print ""
    _fit_ff(anglist, args.labels, ffmol.angles)
