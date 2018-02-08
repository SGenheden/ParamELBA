# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to make a CG structure of an atomistic one.

The force field and mapping needs to be specified in an XML-file.

The output file will be in .gro format.

Examples
--------
    tocg.py -i initial.gro -o cg.gro -x elba.xml
    tocg.py -i initial.gro -o cg.gro -x martini.xml --group 4
"""

import argparse
import sys

import MDAnalysis as md
import MDAnalysis.lib.mdamath as mdutil
import numpy as np

import elbalib

def _reduce_residues(residues, ffmols, size):
    """
    Reduce the number of residue by grouping them together
    """
    #taken = [len(ffmol.beadnames) > 1 for ffmol in ffmols]
    taken = [ffmol.resname != "SOL" for ffmol in ffmols]
    skip = [False for ffmol in ffmols]
    norig = len(residues) - sum(taken)

    ntaken = 0
    i = 0
    while ntaken < norig :
        first_res = residues[taken.index(False)]
        taken[first_res.resid-1] = True
        coord = first_res.atoms.positions
        # Find size - 1 other residues closes to the first non-taken residue
        resid = [res.resid-1 for res, itaken in zip(residues,taken) if not itaken]
        dist = [mdutil.norm(coord[0, :] - res.atoms.positions[0, :])
                for res, itaken in zip(residues,taken) if not itaken]
        sortlst = np.asarray(dist).argsort()
        # Marken the closest residues for skipping
        for isort in sortlst[:size-1] :
            idsort = resid[isort]
            skip[idsort] = True
            taken[idsort] = True
            coord = coord + residues[idsort].atoms.positions
        # Move the first non-taken residue to the centroid of the group
        coord = coord / float(size)
        first_res.atoms.positions = coord
        ntaken += size

    return skip

if __name__ == '__main__':

    print " ".join(sys.argv)

    parser = argparse.ArgumentParser(description="Transform all-atom structure to CG")
    parser.add_argument('-i', '--infile', help="the all-atom input file")
    parser.add_argument('-o', '--outfile', help="the CG output file")
    parser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    parser.add_argument('--group', type=int, help="group water molecules")    
    args = parser.parse_args()

    ff = elbalib.Elba()
    ff.load(args.xml)
    universe = md.Universe(args.infile,args.infile)
    allsel = universe.select_atoms("all")

    # First make sure that each residue can be mapped to CG
    ffmols = []
    for residue in allsel.residues:
        ffmol = ff.find_molecule(residue.resname)
        if ffmol is None:
            raise Exception("Could not find residue %s in XML-file"%residue.name)
        else:
            ffmols.append(ffmol)

    skip = []
    if args.group is not None :
        skip = _reduce_residues(allsel.residues, ffmols, args.group)

    nbeads = 0
    lines = []
    for i, (residue, ffmol) in enumerate(zip(allsel.residues,ffmols)):
        if skip and skip[i] : continue
        atomnames = [a.name.lower() for a in residue.atoms]
        transmat = ffmol.transmat(atomnames)
        xyz = np.dot(transmat,residue.atoms.positions)
        xyz = xyz / 10.0
        for ibead,(beadname,pos) in enumerate(zip(ffmol.beadnames,xyz),nbeads):
            bead = ffmol.beads[beadname]
            lines.append("%5d%5s%5s%5d%8.3f%8.3f%8.3f\n"%(residue.resnum,residue.resname,bead.name,ibead+1,pos[0],pos[1],pos[2]))
        nbeads = nbeads + len(ffmol.beadnames)

    with open(args.outfile,"w") as f:
        f.write("CG of %s\n"%args.infile)
        f.write("%8d\n"%nbeads)
        for line in lines:
            f.write(line)
        f.write("%8.3f%8.3f%8.3f\n"%tuple(universe.trajectory.ts.dimensions[:3]/10.0))
