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

def _calc_bond(xyz):
    return mdutil.norm(xyz[0, :] - xyz[1, :])


def _calc_angle(xyz):
    a = xyz[0, :] - xyz[1, :]
    b = xyz[2, :] - xyz[1, :]
    return np.rad2deg(np.arccos(np.dot(a, b) / (mdutil.norm(a) * mdutil.norm(b))))


def _process_res(residue,processor,**kwargs) :
    """
    Called for each residue at every time step of the trajectory
    Print out the bond and angles to the output file
    """
    f = kwargs["fileobj"]
    f.write("%d"%processor.currtime)
    for i, bondsel in enumerate(residue.bonds):
        f.write("%8.3f" % _calc_bond(residue.xyz[bondsel, :]))
    for i, angsel in enumerate(residue.angles):
        f.write("%8.3f" % _calc_angle(residue.xyz[angsel, :]))
    f.write("\n")


if __name__ == '__main__':

    print " ".join(sys.argv)

    # Setup the trajectory processor
    processor = mdlib.AaCgTrajectoryProcessor("Track CG bond and angle distributions through an MD trajectory")
    processor.argparser.add_argument('-o', '--out', help="the output", default="bondang.txt")
    
    processor.setup()
    # Add bond and angle indices to each of the defined residues
    for res in processor.residues :
        res.bonds = [bond.indices(res.ffmol) for bond in res.ffmol.bonds]
        res.angles = [angle.indices(res.ffmol) for angle in res.ffmol.angles]

    # Open the output file in a context manager and then process the trajectory
    with open(processor.args.out, "w") as f:
        processor.process(None,residuefunc=_process_res,fileobj=f)
                
