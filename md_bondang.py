# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to track CG bond and angles during an MD trajectory.

The MD trajectory can be either CG or AA and in the latter case
a transformation is carried out from the AA coordinates to the CG
representation.

The force field definition, i.e. the bond and angles to be tracked
are stored in an XML format and read by the classes in the elbalib
module.

A reference structure need to be given and this is used to find CG
molecules defined in the XML file. If no additional structure file
is given it is assumed that the trajectory is a CG trajectory and
this is used to read the MD trajectory.

If an additional structure is given it is assumed that the trajectory
is AA and this structure is used to read the MD trajectory as well
as to setup the transformation from the AA to the CG representation.

Examples
--------
md_bondang.py -f prod.dcd -r popc_elba.pdb -x elba.xml
md_bondang.py -f md2_whole.xtc -s md1.gro -r popc_elba.pdb -x elba.xml
"""

import sys

import MDAnalysis.core.util as mdutil
import numpy as np

import mdlib

class BondAngAnalysis(mdlib.AaCgAction):
    """
    Class to analyse bond and angles during a trajectory

    Attributes:
    -----------
    f : FileObject
        the output file
    """
    def add_arguments(self, parser):
        parser.add_argument('-o', '--out', help="the output", default="bondang.txt")

    def setup(self, args):
        self.out = args.out
        self.f = open(self.out, "w")

    def setup_residues(self):
        for res in self.processor.residues :
            res.bonds = [bond.indices(res.ffmol) for bond in res.ffmol.bonds]
            res.angles = [angle.indices(res.ffmol) for angle in res.ffmol.angles]

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        try:
            self.f.close()
        except:
            pass
        return False

    def resprocess(self, residue):
        """
        Calculates all bond and angle values and print them to disc
        """
        self.f.write("%d"%self.processor.currtime)
        for i, bondsel in enumerate(residue.bonds):
            self.f.write("%8.3f" % self._calc_bond(residue.xyz[bondsel, :]))
        for i, angsel in enumerate(residue.angles):
            self.f.write("%8.3f" % self._calc_angle(residue.xyz[angsel, :]))
        self.f.write("\n")

    def _calc_bond(self,xyz):
        return mdutil.norm(xyz[0, :] - xyz[1, :])

    def _calc_angle(self,xyz):
        a = xyz[0, :] - xyz[1, :]
        b = xyz[2, :] - xyz[1, :]
        return np.rad2deg(np.arccos(np.dot(a, b) / (mdutil.norm(a) * mdutil.norm(b))))

if __name__ == '__main__':

    processor = mdlib.AaCgTrajectoryProcessor("Track CG bond and angle distributions through an MD trajectory")
    with BondAngAnalysis(processor) as analysis:
        processor.setup(printargs=True)
        processor.process()
