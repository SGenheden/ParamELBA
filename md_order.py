# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Program to calculate CG tail order parameter from MD trajectory

The MD trajectory can be either CG or AA and in the latter case
a transformation is carried out from the AA coordinates to the CG
representation.

The force field definition are stored in an XML format and read 
by the classes in the elbalib module.

A reference structure need to be given and this is used to find CG
molecules defined in the XML file. If no additional structure file
is given it is assumed that the trajectory is a CG trajectory and
this is used to read the MD trajectory.

If an additional structure is given it is assumed that the trajectory
is AA and this structure is used to read the MD trajectory as well
as to setup the transformation from the AA to the CG representation.

The carbon beads that form the tail on which the order parameters
should be calculated are given as command-line arguments. One or
more chains can be given.

At the moment all the lipids have to be of a single type.

Examples
--------
md_order.py -f md2_whole.xtc -r popc_martini.gro -x martini.xml -c C1A-C2A-C3A-C4A C1B-C2B-D3B-C4B-C5B
md_order.py -f md2_whole.xtc -s md1.gro -r popc_martini.pdb -x martini.xml -c C1A-C2A-C3A-C4A C1B-C2B-D3B-C4B-C5B
"""

import sys

import MDAnalysis.core.util as mdutil
import numpy as np

import mdlib

def _process_res(residue,processor,**kwargs):
    """
    Called for each residue at every time step of the trajectory
    Calculates order parameters for each bond vector in each residue
    """
    for i,chain in enumerate(residue.chains) :
        # Atom2 - Atom1
        vec = residue.xyz[chain[1:],:] - residue.xyz[chain[:-1],:]
        # Projection with normal
        proj = np.multiply(vec,kwargs["norm"]).sum(axis=1)**2 / np.sum(vec**2,axis=1)
        # Order param, these will be printed out by _process_post
        residue.order[i] = 0.5*(3.0*proj-1)


def _process_post(processor,**kwargs):
    """
    Called at every time step of the trajectory, after the residues has been updated
    Calculated mol-average order parameters and print them to the outputfile
    """           
    f = kwargs["fileobj"] 
    f.write("%d"%processor.currtime)
    # Calculates molecular-averaged order parameters
    avorder = [np.zeros(o.shape[0]) for o in processor.residues[0].order]
    for residue in processor.residues :
        for ichain,rorder in enumerate(residue.order) :
            avorder[ichain] += rorder
    # Print out the average to the output file
    one_over_n = 1 / float(len(processor.residues))
    for i, o in enumerate(avorder):
        omolav = o * one_over_n
        f.write(" %s" % " ".join("%.3f"%oo for oo in omolav))
    f.write("\n")


if __name__ == '__main__':

    print " ".join(sys.argv)

    # Setup the trajectory processor
    processor = mdlib.AaCgTrajectoryProcessor("Calculate CG tail order parameter from MD trajectory")
    processor.argparser.add_argument('-c', '--chains',nargs="+", help="the carbon chains")
    processor.argparser.add_argument('-o', '--out', help="the output", default="order.txt")
    
    processor.setup()
    # Add bead indices for each chain to track
    for res in processor.residues :
        res.chains = []
        res.order = []
        for chain in processor.args.chains :
            indices = np.array([res.ffmol.beads[b].idx for b in chain.split("-")],dtype=int)
            res.chains.append(indices)
            res.order.append(np.zeros(indices.shape[0]-1))

    # Assumes that the normal is along the z-axis
    bilayer_norm = np.array([0.0,0.0,1.0])

    # Open the output file in a context manager and then process the trajectory
    with open(processor.args.out, "w") as f:
        processor.process(None,residuefunc=_process_res,postfunc=_process_post,fileobj=f,norm=bilayer_norm)
                
