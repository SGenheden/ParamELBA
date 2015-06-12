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

A current restriction is that the molecules to be tracked is of a
single kind, e.g. a single lipid type in a membrane.

Examples
--------
md_bondang.py -f prod.dcd -r popc_elba.pdb -x elba.xml
md_bondang.py -f md2_whole.xtc -s md1.gro -r popc_elba.pdb -x elba.xml
"""

import argparse
import sys

import MDAnalysis as md
import MDAnalysis.core.util as mdutil
import numpy as np

import elbalib


def _calc_bond(xyz):
    return mdutil.norm(xyz[0, :] - xyz[1, :])


def _calc_angle(xyz):
    a = xyz[0, :] - xyz[1, :]
    b = xyz[2, :] - xyz[1, :]
    return np.rad2deg(np.arccos(np.dot(a, b) / (mdutil.norm(a) * mdutil.norm(b))))


if __name__ == '__main__':

    print " ".join(sys.argv)

    # Command-line input
    parser = argparse.ArgumentParser(description="Track CG bond and angle distributions through an MD trajectory")
    parser.add_argument('-f', '--file', nargs="+", help="the trajectory file.", default="")
    parser.add_argument('-r', '--ref', help="a CG reference file", default="ref.pdb")
    parser.add_argument('-s', '--struct', help="a AA structure file")
    parser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    parser.add_argument('-o', '--out', help="the output", default="bondang.txt")
    parser.add_argument('-n', '--naming', nargs=2, help="naming convention change")
    parser.add_argument('--skip', type=int, help="skip this many snapshots", default=0)
    parser.add_argument('--dt', type=float, help="the number of ps for each snapshot", default=10)
    args = parser.parse_args()

    # Load the force field definitions
    ff = elbalib.Elba()
    ff.load(args.xml)

    # Load a universe for the MD trajectory
    if args.struct is None:
        universe = md.Universe(args.ref, args.file)
    else:
        universe = md.Universe(args.struct, args.file)
    # Load the reference CG universe
    refuni = md.Universe(args.ref)

    # Check if we need to translate the atom names
    namtrans = None
    if args.naming is not None:
        namref = [line.strip() for line in open(args.naming[0], 'r').readlines()]
        nammob = [line.strip() for line in open(args.naming[1], 'r').readlines()]
        namtrans = {}
        for nr, nm in zip(namref, nammob):
            if nm != "*":
                namtrans[nm] = nr

    # Pre-processing
    resids = []
    mdresidues = []
    bond_dict = {}
    ang_dict = {}
    trans_mats = {}
    ffmols = {}
    for r in refuni.selectAtoms("all").residues:
        ffmol = ff.find_molecule(r.name)
        if ffmol is None or not (ffmol.bonds and ffmol.angles):
            continue
        resid = "resid %d" % r.id
        resids.append(resid)
        mdresidues.append(universe.selectAtoms(resid))
        ffmols[resid] = ffmol

        # Setup bond and angle selections that are to be tracked
        bond_dict[resid] = [bond.indices(ffmol) for bond in ffmol.bonds]
        ang_dict[resid] = [angle.indices(ffmol) for angle in ffmol.angles]

        # Setup transformation matrices if necessary
        if args.struct is not None:
            if namtrans is None:
                atomnames = [a.name.lower() for a in mdresidues[-1]]
            else:
                atomnames = [namtrans[a.name].lower() for a in mdresidues[-1]]
            trans_mats[resid] = ffmol.transmat(atomnames)

    # Open the output file in a context manager
    with open(args.out, "w") as f:

        # Loop over the coordinates
        for ti, ts in enumerate(universe.trajectory):

            if ti * args.dt % 1000 == 0:
                print "%d" % (ti * args.dt),
                sys.stdout.flush()
            if ti * args.dt < args.skip:
                continue

            # Loop over all residues
            for ri, (resid, mdres) in enumerate(zip(resids, mdresidues)):
                # Get the coordinates, either out-of-the-box CG or transformed AA
                if args.struct is None:
                    xyz = mdres.get_positions()
                else:
                    xyz = np.dot(trans_mats[resid], mdres.get_positions())
                # Print out bond and angle values
                f.write("%d" % (ti * args.dt))
                for i, bondsel in enumerate(bond_dict[resid]):
                    f.write("%8.3f" % _calc_bond(xyz[bondsel, :]))
                for i, angsel in enumerate(ang_dict[resid]):
                    f.write("%8.3f" % _calc_angle(xyz[angsel, :]))
                f.write("\n")
