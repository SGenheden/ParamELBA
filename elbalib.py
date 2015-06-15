# Author: Samuel Genheden samuel.genheden@gmail.com

"""
Classes to store the definition on an ELBA force field.
It contains classes to store definition of beads, atom types, connectivity and molecules.
These definitions are typically read from an XML file.
"""

import xml.etree.ElementTree as ET

import numpy as np
import scipy.stats as stats


class FfType(object):
    """
    Class to store the definition of an force field type

    Attributes
    ----------
    id : integer
      the identifier of the type
    props : dictionary of floats
      the properties of this type
    """

    def __init__(self):
        self.id = -1  # Unassigned
        self.props = {}

    def __str__(self):
        return "Type %d %s" % (self.id, " ".join("%s=%.3f" % (prop, val) for prop, val in self.props.iteritems()))

    def set(self, **kwargs):
        """
        Sets the property of this type from a dictionary, typically an XML node
        """
        if "id" in kwargs:
            self.id = int(kwargs["id"])
        else:
            return
        for prop in kwargs:
            if prop == "id":
                continue  # Already took care of id above
            self.props[prop] = float(kwargs[prop])
            setattr(self, prop, float(kwargs[prop]))

    def potential(self, coord):
        """
        Return the potential given some coordinate
        Needs to be implemented by sub classes
        """
        return 0.0

    def distribution(self, coord, RT):
        """
        Return the distribution of the potential
        Needs to be implemented by sub classes
        """
        return 0.0

    def statmoments(self, RT):
        """
        Return the standard deviation and the mean
        Needs to be implemented by sub classes
        """
        return 0.0, 0.0


class BondType(FfType):
    """
    Class to store the definition of a bond type
    """

    def potential(self, coord):
        return self.k * self.kf * (coord - self.r0) * (coord - self.r0)

    def distribution(self, coord, RT):
        mean, std = self.statmoments(RT)
        return stats.norm.pdf(coord, loc=mean, scale=std)

    def statmoments(self, RT):
        return self.r0, np.sqrt(RT / (2.0 * self.k * self.kf))


class AngleType(FfType):
    """
    Class to store the definition of an angle type
    """

    def potential(self, coord):
        costheta0 = np.cos(np.deg2rad(self.theta0))
        coscoord = np.cos(np.deg2rad(coord))
        return self.k * self.kf * (coscoord - costheta0) * (coscoord - costheta0)

    def distribution(self, coord, RT):
        if coord.min() < -1 or coord.max() > 1.0:
            coscoord = np.cos(np.deg2rad(coord))
        else:
            coscoord = coord
        mean, std = self.statmoments(RT)
        return stats.norm.pdf(coscoord, loc=mean, scale=std)

    def statmoments(self, RT):
        costheta0 = np.cos(np.deg2rad(self.theta0))
        return costheta0, np.sqrt(RT / (2.0 * self.k * self.kf))


class ElbaBead(object):
    """
    Class to store the definition of an ELBA bead

    Attributes
    ----------
    name : string
      the identifier of the bead
    beadtype : BeadType
      the type definition
    idx : integer
      the index of the bead in the molecule that contains it
    mapping : list of string
      the all-atom mapping of this bead
    """

    def __init__(self):
        self.name = None  # Unassigned
        self.beadtype = None
        self.idx = -1
        self.mapping = []

    def __str__(self):
        if self.name is not None:
            if isinstance(self.beadtype, FfType):
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


class Connectivity(object):
    """
    Class to store molecular connectivity

    Attributes
    ----------
    contype : FfType
      the connectivity force field type
    beads : list of ElbaBead
      the beads in the connectivity
    """

    def __init__(self):
        self.contype = -1  # Unassigned
        self.beads = []

    def __str__(self):
        typestr = "%d" % self.type if isinstance(self.type, int) else "(" + self.type.__str__() + ")"
        return "Connectivity type=%s, %s" % (typestr, self.atomstring())

    def atomstring(self):
        return "-".join(b.name if isinstance(b, ElbaBead) else b for b in self.beads)

    def parse(self, element, molecule, typelist):
        """
        Parse connectivity from an XML tree

        Attributes
        ----------
        element : xml.etree.Element
          the node to parse on the XML tree
        molecule : ElbaMolecule 
          the molecule that has this connectivity
        typelist : dictionary of FfType
          the list of types for this connectivity
        """
        if "type" in element.attrib:
            self.type = int(element.attrib["type"])
            if self.type in typelist:
                self.type = typelist[self.type]
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
                if child.tag != "bead":
                    continue
                b = ElbaBead()
                b.parse(child, parent)
                if b.name is not None:
                    b.idx = len(self.beadnames)
                    self.beadnames.append(b.name)
                    self.beads[b.name] = b

        # Then find the connectivity
        tag2list = {"bond": self.bonds, "angle": self.angles, "dipolerest": self.dipolerestraints}
        parentlist = {"bond": parent.bondtypes, "angle": parent.angletypes, "dipolerest": parent.dipoletypes}
        for conroot in element.findall("./connectivity"):
            for child in conroot:
                if child.tag in tag2list:
                    tag2list[child.tag].append(Connectivity())
                    tag2list[child.tag][-1].parse(child, self, parentlist[child.tag])

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
                self._transmat[i, idx] = 1 / float(counts[idx]) * one_over_n  # This does not work!
                self._transmat[i, idx] = one_over_n
        return self._transmat

    def write_cgtools_input(self):

        print "[maptype]"
        print "GC\n"

        print "[residues]"
        print self.resname + "\n"

        print "[frames]"
        print "-1\n"

        print "[mapping]"
        for bname in self.beadnames:
            bead = self.beads[bname]
            print "%s B%d %s" % (bname, bead.beadtype.id if isinstance(bead.beadtype, FfType) else bead.beadtype, " ".join([m.upper() for m in bead.mapping]))
        print ""

        print "[length]"
        for bond in self.bonds:
            print "%s %s" % (bond.beads[0].name, bond.beads[1].name)
        print ""

        print "[angle]"
        for angle in self.angles:
            print "%s %s %s" % (angle.beads[0].name, angle.beads[1].name, angle.beads[2].name)
        print ""


class Elba(object):
    """
    Class to store the definition of the ELBA force field

    Attributes
    ----------
    molecules : dictionary of ElbaMolecule
      the defined molecules
    beadtypes : dictionary of FfType
      the defined beadtypes
    bondtypes : dictionary of BondType
      the defined bondtypes
    angletypes : dictionary of AngleType
      the defined angletypes
    dipoletypes : dictionary of FfType
      the defined dipole restraint types
    units : string
      identifier for the units
    """

    def __init__(self):
        self.molecules = {}
        self.beadtypes = {}
        self.bondtypes = {}
        self.angletypes = {}
        self.dipoletypes = {}

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

        # Parse bead, bond angle angle types
        typedict = {"beadtype": FfType, "bondtype": BondType, "angletype": AngleType, "dipoletype": FfType}
        listdict = {"beadtype": self.beadtypes, "bondtype": self.bondtypes,
                    "angletype": self.angletypes, "dipoletype": self.dipoletypes}
        for child in tree.getroot():
            if child.tag in typedict:
                t = typedict[child.tag]()
                t.set(**child.attrib)
                if t.id > -1:
                    listdict[child.tag][t.id] = t

        # Parse molecules
        for child in tree.getroot():
            if child.tag != "mol":
                continue
            mol = ElbaMolecule()
            mol.parse(child, self)
            if mol.resname is not None:
                self.molecules[mol.resname] = mol


def test(filename):
    elbaobj = Elba()
    elbaobj.load(filename)
    for id, bt in elbaobj.beadtypes.iteritems():
        print bt
    for id, bt in elbaobj.bondtypes.iteritems():
        print bt
        print bt.potential(4.0)
    for id, bt in elbaobj.angletypes.iteritems():
        print bt
        print bt.potential(120.0)

    for id, mol in elbaobj.molecules.iteritems():
        print mol
        for b in mol.bonds:
            print b
        for a in mol.angles:
            print a
        for d in mol.dipolerestraints:
            print d

    print elbaobj.molecules["POPC"].write_cgtools_input()

if __name__ == '__main__':
    import sys
    test(sys.argv[1])
