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

class TailOrderAnalysis(mdlib.AaCgAction):
    """
    Class to analyse tail order parameters during a trajectory

    Attributes:
    -----------
    carbonsel : MDAnalysis.AtomGroup
        the carbon selections
    hydrogensel : MDAnalysis.AtomGroup
        the hydrogen selections
    f_cc : FileObject
        the output file of C-C order params
    f_ch : FileObject
        the output file of C-H order params
    """
    def add_arguments(self, parser):
        parser.add_argument('-c', '--carbons',nargs="+", help="the carbon chains")
        parser.add_argument('-y', '--hydrogens',nargs="+", help="the hydrogen names")
        parser.add_argument('-x', '--extranames',nargs="+", help="extra hydrogen names for non-serial names")
        parser.add_argument('-o', '--out', help="the output prefix", default="order_aa")

    def setup(self, args):
        self.carbonsel = []
        self.hydrogensel = []
        extranames = args.extranames[::-1]
        for carbon,hydrogen in  zip(args.carbons,args.hydrogens):
            atoms = carbon.split("-")
            if len(atoms) == 1 : atoms = self._expand_selection(carbon)
            atomsels = [self.processor.universe.selectAtoms("name %s"%atom) for atom in atoms]
            self.carbonsel.append(atomsels)

            atoms = hydrogen.split("-")
            atomsels = []
            if len(atoms) == 1 : atoms = self._expand_selection(hydrogen)
            for atom in atoms:
                sel = self.processor.universe.selectAtoms("name %s"%atom)
                if len(sel) == 0:
                    sel = self.processor.universe.selectAtoms("name %s"%extranames.pop())
                atomsels.append(sel)
            self.hydrogensel.append(atomsels)

        self.out = args.out
        # Assumes that the normal is along the z-axis
        self.normal = np.array([0.0,0.0,1.0])

    def __enter__(self):
        self.f_cc = open(self.out+"_cc.txt", "w")
        self.f_ch = open(self.out+"_ch.txt", "w")
        return self

    def __exit__(self, type, value, traceback):
        self.f_cc.close()
        self.f_ch.close()
        return False

    def process(self):
        self.f_cc.write("%d"%self.processor.currtime)
        self.f_ch.write("%d"%self.processor.currtime)
        for carbonsels,hydrogensels in zip(self.carbonsel,self.hydrogensel) :
            for a1,a2 in zip(carbonsels[:-1],carbonsels[1:]):
                self.f_cc.write(" %.3f"%self._calc_order(a1,a2,self.normal))
            for a1,a2 in zip(carbonsels[:-1],hydrogensels[:-1]):
                self.f_ch.write(" %.3f"%self._calc_order(a1,a2,self.normal))
        self.f_cc.write("\n")
        self.f_ch.write("\n")

    def _calc_order(self,a1,a2,norm):
        # Atom2 - Atom1
        vec = a2.get_positions() - a1.get_positions()
        # Projection with normal
        proj = np.multiply(vec,norm).sum(axis=1)**2 / np.sum(vec**2,axis=1)
        # Order param
        return np.abs(0.5*(3.0*proj.mean()-1))

    def _expand_selection(self,selection):
        """
        Expand a simple pattern with {A..B} format to an integer range
        """
        first = selection.index("{")
        last  = selection.index("}")
        start,end = map(int,selection[first+1:last].split(".."))
        return ["%s%d%s"%(selection[:first],i,selection[last+1:]) for i in range(start,end+1)]

if __name__ == '__main__':

    processor = mdlib.TrajectoryProcessor("Calculate tail order parameter from MD trajectory")
    with TailOrderAnalysis(processor) as analysis:
        processor.setup(printargs=True)
        processor.process()
