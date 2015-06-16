# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to help with the processing of MD trajectories
"""

import argparse
import sys

import MDAnalysis as md
import numpy as np

import elbalib


class TrajectoryProcessor(object):
    """
    Class to process an MD trajectory
    
    The program that uses this initialises an instance that setups up an
    argparse command-line interpreter.
    The program can add its own arguments to this parser before calling setup()
    that parse the arguments and initialises the md universe object
    To analyse the trajectory the program then calls the process() routine and
    passes a function that will be called for each snapshot in the trajectory
    
    Attributes
    ----------
    argparser : argparse.ArgumentParser
        the command line argument parser
    args : argparse.Namespace
        the parsed arguments
    currtime : float
        the current time in ps when processing the trajectory
    namtrans : dictionary of strings
        name translations
    skip : float
        the number of ps to skip
    dt : float
        the number of ps per snapshot
    universe : MDAnalysis.Universe
        the MDAnalysis universe holding the trajectory and its state
    """
    
    
    def __init__(self, description):
        self.argparser = argparse.ArgumentParser(description=description)
        self.argparser.add_argument('-f', '--file', nargs="+", help="the trajectory file.")
        self.argparser.add_argument('-s', '--struct', help="a structure file")
        self.argparser.add_argument('-n', '--naming', nargs=2, help="naming convention change") 
        self.argparser.add_argument('--skip', type=int, help="skip this many snapshots", default=0)
        self.argparser.add_argument('--dt', type=float, help="the number of ps for each snapshot", default=10)
        self.universe = None
        self.namtrans = None
        self.currtime = 0
        self.dt = 0
        self.skip = 0
        
    def setup(self):
        """
        Parsing command-line arguments and setting up the MD universe
        This routine should be called before process()
        """
        self.args = self.argparser.parse_args()
        self.skip = self.args.skip
        self.dt = self.args.dt
        self.setup_universe()

        # Check if we need to translate the atom names
        self.namtrans = None
        if self.args.naming is not None:
            namref = [line.strip() for line in open(self.args.naming[0], 'r').readlines()]
            nammob = [line.strip() for line in open(self.args.naming[1], 'r').readlines()]
            self.namtrans = {}
            for nr, nm in zip(namref, nammob):
                if nm != "*":
                    self.namtrans[nm] = nr  
    
    def setup_universe(self):
        if self.args.file is None or self.args.struct is None :
            raise Exception("Both a universe and structure file needs to be specified")
        self.universe = md.Universe(self.args.struct, self.args.file)
        
    def process(self,processfunc,**kwargs):
        """
        Loop over the MD trajectory and process it
        
        The processfunc will be called at every step of the trajectory and
        it should take a TrajectoryProcessor instance (self) as well and a
        dictionary of additional arguments.
        
        Arguments
        ---------
        processfunc : function handle
            a function that will be called at each snapshot
        **kwargs : dictionary
            a number of extra arguments passed to the process function
        """
        if self.universe is None : return
        
        # Loop over the coordinates
        for ti, ts in enumerate(self.universe.trajectory):
            self.currtime = ti*self.dt
            if self.currtime % 1000 == 0:
                print "%d" % (self.currtime),
                sys.stdout.flush()
            if self.currtime < self.skip:
                continue
                
            processfunc(self,**kwargs)


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


class AaCgTrajectoryProcessor(TrajectoryProcessor):
    """
    A processor for a trajectory that can be either AA or CG
    
    Attributes
    ----------
    ff : elbalib.Elba
        the force field definitions
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
        self.ff = None
        self.refuniverse = None
        self.residues = []
        self.iscgtraj = False
        
    def setup(self):
        super(AaCgTrajectoryProcessor,self).setup()
        
        # Load the force field definitions
        self.ff = elbalib.Elba()
        self.ff.load(self.args.xml)
        
        # Load the reference CG universe
        self.refuniverse = md.Universe(self.args.ref)  
        
        # Setup of a list of residues defined in the force field that can be processed           
        self.residues = []
        for r in self.refuniverse.selectAtoms("all").residues:
            ffmol = self.ff.find_molecule(r.name)
            if ffmol is None or not (ffmol.bonds and ffmol.angles):
                continue
            self.residues.append(_Residue(r.id, ffmol, self.universe))
            if not self.iscgtraj : self.residues[-1].make_transmat(self.namtrans)
                   
    def setup_universe(self):
        self.iscgtraj = self.args.struct is None
        if self.iscgtraj:
            self.args.struct = self.args.ref
        super(AaCgTrajectoryProcessor,self).setup_universe()
        
    def process(self,processfunc,**kwargs):
        """
        Loop over the md trajectory and update the coordinates of the stored residues

        The processfunc argument is ignored and the program using this class should
        pass a function to one or more of prefunc, residuefunc or postfunc to perform
        the actual processing. prefunc is called before the coordinates of the
        residues are updated, residuefunc is called once for each residue and
        postfunc is called after all residues has been updated.

        The other arguments in **kwargs are passed along to these functions
        """

        def _total_process(processor,**kwargs):
            if "prefunc" in kwargs and kwargs["prefunc"] is not None:
                kwargs["prefunc"](processor,**kwargs)
            # Loop over all residues
            for res in self.residues:
                res.update_coordinates()
                if "residuefunc" in kwargs and kwargs["residuefunc"] is not None:
                    kwargs["residuefunc"](res,processor,**kwargs)
            if "postfunc" in kwargs and kwargs["postfunc"] is not None:
                kwargs["postfunc"](processor,**kwargs)
        
        super(AaCgTrajectoryProcessor,self).process(_total_process,**kwargs)
