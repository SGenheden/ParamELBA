# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate CG tail order parameter from MD trajectory

It will calculate both C-C and C-H order parameters

The carbon atoms that form the tail on which the order parameters
should be calculated are given as command-line arguments. One or
more chains can be given.

The hydrogen atoms are also given as a command-line argument.

At the moment all the lipids have to be of a single type.

Examples
--------
md_order_aa.py -f md2_whole.xtc -s md1.gro -c C2{2..18} -y H{2..18}S -x H91 H101
"""

import sys

import MDAnalysis.core.util as mdutil
import numpy as np

import mdlib

def _doprocess(processor,**kwargs):
    """
    Called at every time step of the trajectory
    Calculates order parameters for each bond vector and prints it to output
    """

    def _calc_order(a1,a2,norm):
        # Atom2 - Atom1
        vec = a2.get_positions() - a1.get_positions()
        # Projection with normal
        proj = np.multiply(vec,norm).sum(axis=1)**2 / np.sum(vec**2,axis=1)
        # Order param
        return np.abs(0.5*(3.0*proj.mean()-1))

    f_cc = kwargs["fileobj_cc"] 
    f_ch = kwargs["fileobj_ch"] 
    f_cc.write("%d"%processor.currtime)
    f_ch.write("%d"%processor.currtime)
    for carbonsels,hydrogensels in zip(processor.carbonsel,processor.hydrogensel) :
        for a1,a2 in zip(carbonsels[:-1],carbonsels[1:]):
            f_cc.write(" %.3f"%_calc_order(a1,a2,kwargs["norm"]))
        for a1,a2 in zip(carbonsels[:-1],hydrogensels[:-1]):
            f_ch.write(" %.3f"%_calc_order(a1,a2,kwargs["norm"]))
    f_cc.write("\n")
    f_ch.write("\n")
    

def _expand_selection(selection):
    """
    Expand a simple pattern with {A..B} format to an integer range
    """
    first = selection.index("{")
    last  = selection.index("}")
    start,end = map(int,selection[first+1:last].split(".."))
    return ["%s%d%s"%(selection[:first],i,selection[last+1:]) for i in range(start,end+1)]


if __name__ == '__main__':

    print " ".join(sys.argv)

    # Setup the trajectory processor
    processor = mdlib.TrajectoryProcessor("Calculate tail order parameter from MD trajectory")
    processor.argparser.add_argument('-c', '--carbons',nargs="+", help="the carbon chains")
    processor.argparser.add_argument('-y', '--hydrogens',nargs="+", help="the hydrogen names")
    processor.argparser.add_argument('-x', '--extranames',nargs="+", help="extra hydrogen names for non-serial names")
    processor.argparser.add_argument('-o', '--out', help="the output prefix", default="order_aa")
    
    processor.setup()
    # Add atom selections to the processor
    processor.carbonsel = []
    processor.hydrogensel = []
    extranames = processor.args.extranames[::-1]
    for carbon,hydrogen in  zip(processor.args.carbons,processor.args.hydrogens):
        atoms = carbon.split("-")
        if len(atoms) == 1 : atoms = _expand_selection(carbon)
        atomsels = [processor.universe.selectAtoms("name %s"%atom) for atom in atoms]
        processor.carbonsel.append(atomsels)

        atoms = hydrogen.split("-")
        atomsels = []
        if len(atoms) == 1 : atoms = _expand_selection(hydrogen)
        for atom in atoms:
            sel = processor.universe.selectAtoms("name %s"%atom)
            if len(sel) == 0:
                sel = processor.universe.selectAtoms("name %s"%extranames.pop())    
            atomsels.append(sel)
        processor.hydrogensel.append(atomsels)
        

    # Assumes that the normal is along the z-axis
    bilayer_norm = np.array([0.0,0.0,1.0])

    # Open the output files in a context manager and then process the trajectory
    with open(processor.args.out+"_cc.txt", "w") as f_cc, open(processor.args.out+"_ch.txt", "w") as f_ch :
        processor.process(_doprocess,fileobj_cc=f_cc,fileobj_ch=f_ch,norm=bilayer_norm)
                
