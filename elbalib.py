# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to store the definition on an ELBA force field.
It contains classes to store definition of beads, atom types, connectivity and molecules.
These definitions are typically read from an XML file.
"""

import xml.etree.ElementTree as ET

import numpy as np


class ElbaBead(object):
    """
    Class to store the definition of an ELBA bead

    Attributes
    ----------
    name : string
      the identifier of the bead
    beadtype : BeadType
      the type definition
    mapping : list of string
      the all-atom mapping of this bead
    """

    def __init__(self):
        self.name = None  # Unassigned
        self.beadtype = None
        self.mapping = []

    def __str__(self):
        if self.name is not None:
            if isinstance(self.beadtype, BeadType):
                return "Bead %s type=(%s) mapping=(%s)" % (self.name, self.beadtype.__str__(), " ".join(self.mapping))
            else:
                return "Bead %s type=%d mapping=(%s)" % (self.name, self.beadtype, " ".join(self.mapping))
        else:
            return "Unassigned bead"

    def parse(self, element, parent):
        """
        Parse bead definition from an XML tree

        Attributes
        ----------
        element : xml.etree.Element
          the node to parse on the XML tree
        parent : Elba
          the owner of this bead
        """
        if "name" in element.attrib:
            self.name = element.attrib["name"]
        else:
            return
        if "type" in element.attrib:
            t = int(element.attrib["type"])
            if t in parent.beadtypes:
                self.beadtype = parent.beadtypes[t]
            else:
                self.beadtype = t
        for m in element.findall("./mapping"):
            self.mapping.extend(m.text.strip().split())


class BeadType(object):
    """
    Class to store the definition of an ELBA bead type

    Attributes
    ----------
    id : integer
      the identifier of the type
    charge : float
      the charge of the bead
    dipole : float
      the magnitude of the dipole vector
    density : float
      the density of the bead
    radius : float
      the radius of the bead
    """

    def __init__(self):
        self.id = -1  # Unassigned
        self.charge = 0.0
        self.dipole = 0.0
        self.density = 0.0
        self.radius = 0.0

    def __str__(self):
        return "Type %d charge=%.2f dipole=%.2f density=%.2f radius=%.2f" % (self.id, self.charge, self.dipole, self.density, self.radius)

    def set(self, **kwargs):
        """
        Sets the property of this bead type from a dictionary, typically an XML node
        """
        if "id" in kwargs:
            self.id = int(kwargs["id"])
        else:
            return
        for prop in ["charge", "dipole", "density", "radius"]:
            if prop in kwargs:
                setattr(self, prop, float(kwargs[prop]))


class Connectivity(object):
    """
    Class to store molecular connectivity

    Attributes
    ----------
    contype : integer
      the connectivity force field type
    beads : list of ElbaBead
      the beads in the connectivity
    """

    def __init__(self):
        self.contype = -1  # Unassigned
        self.beads = []

    def __str__(self):
        return "-".join(b.name if isinstance(b, ElbaBead) else b for b in self.beads)

    def parse(self, element, molecule, parent):
        """
        Parse connectivity from an XML tree

        Attributes
        ----------
        element : xml.etree.Element
          the node to parse on the XML tree
        molecule : ElbaMolecule 
          the molecule that has this connectivity
        parent : Elba
          the owner of this molecule
        """
        if "type" in element.attrib:
            self.type = int(element.attrib["type"])
            for atom in element.text.strip().split():
                if atom in molecule.beadnames:
                    self.beads.append(molecule.beads[atom])
                else:
                    self.beads.append(atom)

    def selection(self):
        return " or ".join("name %s" % b.name if isinstance(b, ElbaBead) else "name %s" % b for b in self.beads)

    def indices(self, molecule):
        return np.array([molecule.beadnames.index(b.name) if isinstance(b, ElbaBead) else molecule.beadnames.index(b) for b in self.beads])


class ElbaMolecule(object):
    """
    Class to store the definition of an ELBA molecule

    Attributes
    ----------
    beadnames : list of string
      the name of all beads
    beads : dictionary of ElbaBead
      the beads in the molecule
    bonds : list of Connectivity
      the bonds in the molecule
    angles : list of Connectivity
      the angles in the molecule
    dipolerestraints : list of Connectivity
      the dipolar restraints in this molecule
    resname : string
      the name of the molecule
    mappingff : string
      identifier for the all-atom force field for mapping
    """

    def __init__(self):
        self.beadnames = []
        self.beads = {}
        self.bonds = []
        self.angles = []
        self.dipolerestraints = []
        self.resname = None  # Unassigned
        self._transmat = None

    def __str__(self):
        return "\n".join([self.beads[bn].__str__() for bn in self.beadnames])

    def parse(self, element, parent):
        """
        Parse molecule from an XML tree

        Attributes
        ----------
        element : xml.etree.Element
          the node to parse on the XML tree
        parent : Elba
          the owner of this molecule
        """
        if "resname" in element.attrib:
            self.resname = element.attrib["resname"]
        else:
            return

        if "mappingff" in element.attrib:
            self.mappingff = element.attrib["mappingff"]

        # First find the beads
        for beadsroot in element.findall("./beads"):
            for child in beadsroot:
                if child.tag != "bead": continue
                b = ElbaBead()
                b.parse(child, parent)
                if b.name is not None:
                    self.beadnames.append(b.name)
                    self.beads[b.name] = b

        # Then find the connectivity
        tag2list = {"bond": self.bonds, "angle": self.angles, "dipolerest": self.dipolerestraints}
        for conroot in element.findall("./connectivity"):
            for child in conroot:
                if child.tag in tag2list:
                    tag2list[child.tag].append(Connectivity())
                    tag2list[child.tag][-1].parse(child, self, parent)

    def transmat(self, atomnames, usechache=True):
        if usechache and self._transmat is not None:
            return self._transmat

        # Calculates the number of times an atom appears in a mapping list
        counts = [0] * len(atomnames)
        for bname, bead in self.beads.iteritems():
            for aname in bead.mapping:
                if aname in atomnames:
                    idx = atomnames.index(aname)
                    counts[idx] = counts[idx] + 1
                else:
                    raise Exception("Found mapping atom (%s) not in atom name list" % aname)

        # Now build the transformation matrix
        self._transmat = np.zeros([len(self.beadnames), len(atomnames)])
        for i, bname in enumerate(self.beadnames):
            bead = self.beads[bname]
            one_over_n = 1 / float(len(bead.mapping))
            for aname in bead.mapping:
                idx = atomnames.index(aname)
                self._transmat[i, idx] = 1 / float(counts[idx]) * one_over_n
        return self._transmat


class Elba(object):
    """
    Class to store the definition of the ELBA force field

    Attributes
    ----------
    molecules : dictionary of ElbaMolecule
      the defined molecules
    beadtypes : dictionary of BeadTypes
      the defined beadtypes
    units : string
      identifier for the units
    """

    def __init__(self):
        self.molecules = {}
        self.beadtypes = {}

    def find_molecule(self, name):
        for mol in self.molecules:
            if mol.lower() == name.lower() or mol[:3].lower() == name.lower():
                return self.molecules[mol]
        return None

    def load(self, filename):
        """
        Open an XML definition of the ELBA force field and parse it into python objects
        """
        tree = ET.parse(filename)

        if "units" in tree.getroot().attrib:
            self.units = tree.getroot().attrib["units"]

        # Parse bead types
        for child in tree.getroot():
            if child.tag != "beadtype": continue
            bt = BeadType()
            bt.set(**child.attrib)
            if bt.id > -1:
                self.beadtypes[bt.id] = bt

        # Parse molecules
        for child in tree.getroot():
            if child.tag != "mol": continue
            mol = ElbaMolecule()
            mol.parse(child, self)
            if mol.resname is not None:
                self.molecules[mol.resname] = mol


def test():
    elbaobj = Elba()
    elbaobj.load('elba.xml')
    for id, bt in elbaobj.beadtypes.iteritems():
        print bt

    for id, mol in elbaobj.molecules.iteritems():
        print mol
        for b in mol.bonds:
            print b
        for a in mol.angles:
            print a
        for d in mol.dipolerestraints:
            print d

    names = [line.strip().split()[0].lower() for line in open("popc.xyz", "r").readlines()[2:]]
    xyz = np.array([line.strip().split()[1:] for line in open("popc.xyz", "r").readlines()[2:]], dtype=float)
    t = elbaobj.molecules["POPC"].transmat(names).T
    for i, col in enumerate(t):
        print "%5s" % names[i],
        for d in col:
            print "%.2f" % d,
        print ""
    xyz2 = np.dot(t.T, xyz)
    print xyz2.shape[0]
    print "0"
    for bname, c in zip(elbaobj.molecules["POPC"].beadnames, xyz2):
        print "%s %.3f %.3f %.3f" % (bname, c[0], c[1], c[2])
