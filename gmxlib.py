# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to help with the processing of Gromacs itp-files
"""

import re

import elbalib


class AtomType(elbalib.FfType):
    """
    Class to store the properties of an atom type in a itp-file
    
    Attributes
    ----------
    name : string
      a string identifier
    atomnum : int
      the atomic number
    mass : float
      the mass 
    sigma : float
      the LJ sigma parameter
    epsilon : float
      the LJ epsilon parameter
    
    """

    def read_itp(self,id,line):
        """
        Read a line from an itp file that defines the atom type
        """
        if line.find(";") > -1:
            line,comment = line.split(";")
        cols = line.strip().split()
        self.set(id=id,name=cols[0],atomnum=cols[1],
                    mass=float(cols[2]),sigma=float(cols[5])*10.0,
                    epsilon=float(cols[6])/4.184)

class ItpAtom(elbalib.FfType):
    """
    Class to store the properties of an atom in an ipt-file

    Attributes
    ----------
    serial : int
      the serial number of the atom in the parent list
    type : AtomType
      the atom type
    name : string
      a string identifier
    charge : float
      the charge of the atom
    mass : float
      the mass of the atom
    """

    def read_itp(self,line):
        """
        Read a line from an itp file that defines the atom 
        """
        if line.find(";") > -1:
            line,comment = line.split(";")
        cols = line.strip().split()
        self.set(id=int(cols[0]),type=cols[1],name=cols[4],
                    charge=float(cols[6]),mass=float(cols[7]))

class SimpleItp(object):
    """
    Class to hold the definition of a single molecule consisting of a
    single residue in a Itp-file

    Attributes
    ----------
    atoms : list of ItpAtom
        the atoms in the itp file
    bonds : list of tuples of IptAtom
        the bonds in the itp file
    name : string
        the residue name
    """
    
    def __init__(self):
        self.atoms = []
        self.bonds = []
        self.name = ""

    def _find_section(self,fileobj,section,line):
        match = re.search("\[ *%s *\]"%section,line)
        while line and not match:
            line = fileobj.readline()
            match = re.search("\[ *%s *\]"%section,line)  
        if not line : raise Exception("Could not find section %s in this ITP file"%section)
        line = fileobj.readline()
        while line and line.strip()[0] == ";":
            line = fileobj.readline() 
        return line  

    def bonded_to(self,atom):

        lst = []
        for bond in self.bonds:
            if atom == bond[0]:
                lst.append(bond[1])
            elif atom == bond[1]:
                lst.append(bond[0])
        return lst

    def find_atom(self,name):

        for a in self.atoms:
            if a.name.lower() == name.lower() :
                return a
        return None

    def read_itp(self,filename):           
        """
        Read an itp file from disc
        """

        with open(filename,"r") as f:

            line = f.readline()
            
            # First find the name of this molecule
            line = self._find_section(f,"moleculetype",line)
            cols = line.strip().split()
            self.name = cols[0]

            # Then atoms
            line = self._find_section(f,"atoms",line)
            while line:
                if len(line) > 5 and line.strip()[0] != ";" : 
                    if line.strip()[0] == "[" : break
                    self.atoms.append(ItpAtom())
                    self.atoms[-1].read_itp(line)
                line = f.readline()

            # Then bonds
            line = self._find_section(f,"bonds",line)
            while line:
                if len(line) > 5 and line.strip()[0] != ";"  : 
                    if line.strip()[0] == "[" : break
                    cols = line.strip().split()
                    atm1,atm2 = map(int,cols[:2])
                    self.bonds.append((self.atoms[atm1-1],self.atoms[atm2-1]))
                line = f.readline()

    def read_atomtypes(self,filename):
        """
        Read atomtypes from an itp file and assign AtomType objects to the atoms
        """

        types = {}
        
        with open(filename,"r") as f:

            line = f.readline()
            line = self._find_section(f,"atomtypes",line)
            while line:
                if len(line) > 5 and line.strip()[0] != ";" : 
                    if line.strip()[0] == "[" : break
                    t = AtomType()
                    t.read_itp(len(types),line)
                    types[t.name] = t
                line = f.readline()    
    
        for a in self.atoms:
            if a.type in types : 
                a.type = types[a.type]



                