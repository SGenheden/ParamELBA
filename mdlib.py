# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to help with the processing of MD trajectories
"""

import argparse

import MDAnalysis as md
import numpy as np

import elbalib
from sgenlib import moldyn
from sgenlib import mdactions

class _Residue(object):
    """
    Class to store residue information, used by the AaCgTrajectoryProcessor class

    Attributes
    ----------
    ffmol : elbalib.ElbaMolecule
        the force field molecule corresponding to this residue
    mdsel : MDAnalysis.AtomSelection
        the MD universe selection corresponding to this residue
    resid : string
        the identifier of this string
    transmat : numpy.ndarray
        transformation from AA to CG coordinates for this residue
    xyz : numpy.ndarray
        the current coordinates of this residue
    """

    def __init__(self, id, ffmol, mduniverse):
        self.resid = "resid %d" % id
        self.mdsel = mduniverse.selectAtoms(self.resid)
        self.ffmol = ffmol
        self.transmat = None
        self.xyz = None

    def make_transmat(self,namtrans):
        """
        Creates and store a transformation matrix from AA to CG for this residue
        """
        if namtrans is None:
            atomnames = [a.name.lower() for a in self.mdsel]
        else:
            atomnames = [namtrans[a.name].lower() for a in self.mdsel]
        self.transmat = self.ffmol.transmat(atomnames)

    def update_coordinates(self):
        """
        Updates the xyz property by using the coordinates from mdsel,
        if necessary transforms the coordinates
        """
        if self.transmat is None:
            self.xyz = self.mdsel.get_positions()
        else:
            self.xyz = np.dot(self.transmat,self.mdsel.get_positions())


class AaCgTrajectoryProcessor(moldyn.TrajectoryProcessor):
    """
    A processor for a trajectory that can be either AA or CG

    Attributes
    ----------
    ff : elbalib.Elba
        the force field definitions
    namtrans : dictionary of strings
        name translations
    iscgtraj : boolean
        if this is CG trajectory
    refuniverse : MDAnalysis.Universe
        the reference CG universe
    residues : list of _Residue
        residues to keep track of
    """


    def __init__(self, description):
        super(AaCgTrajectoryProcessor,self).__init__(description)
        self.argparser.add_argument('-r', '--ref', help="a CG reference file", default="ref.pdb")
        self.argparser.add_argument('-x', '--xml', help="an XML file with force field definitions")
        self.argparser.add_argument('-n', '--naming', nargs=2, help="naming convention change")
        self.argparser.add_argument('-m','--mol',help="the name of the molecule")
        self._argnames.extend(["ref", "xml", "naming","mol"])
        self.namtrans = None
        self.ff = None
        self.refuniverse = None
        self.residues = []
        self.iscgtraj = False

    def setup(self,printargs=False):
        super(AaCgTrajectoryProcessor,self).setup(printargs)

        self.ff = elbalib.Elba()
        self.ff.load(self.args.xml)

        self.refuniverse = md.Universe(self.args.ref)

        # Check if we need to translate the atom names
        self.namtrans = None
        if self.args.naming is not None:
            namref = [line.strip() for line in open(self.args.naming[0], 'r').readlines()]
            nammob = [line.strip() for line in open(self.args.naming[1], 'r').readlines()]
            self.namtrans = {}
            for nr, nm in zip(namref, nammob):
                if nm != "*":
                    self.namtrans[nm] = nr

        # Setup of a list of residues defined in the force field that can be processed
        self.residues = []
        for r in self.refuniverse.selectAtoms("all").residues:
            if self.args.mol is not None and r.name.lower() != self.args.mol.lower():
                continue
            ffmol = self.ff.find_molecule(r.name)
            if ffmol is None or not ffmol.beadnames:
                continue
            self.residues.append(_Residue(r.id, ffmol, self.universe))
            if not self.iscgtraj : self.residues[-1].make_transmat(self.namtrans)

        for action in self.actions:
            action.setup_residues()

    def setup_universe(self):
        self.iscgtraj = self.args.struct is None
        if self.iscgtraj:
            self.args.struct = self.args.ref
        super(AaCgTrajectoryProcessor,self).setup_universe()

    def call_actions(self):
        for action in self.actions:
            action.preprocess()

        for res in self.residues:
            res.update_coordinates()
            for action in self.actions:
                action.resprocess(res)


        for action in self.actions:
            action.process()
            if self.dosubsample and action.dosubsample and \
                    self.currtime % self.freq == 0:
                action.subsample()

class AaCgAction(mdactions.TrajectoryAction):

    def setup_residues():
        pass

    def preprocess(self):
        pass

    def resprocess(self,res):
        pass
