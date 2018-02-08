# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse
import copy
import sys

import numpy as np

import elbalib
import gmxlib


if __name__ == '__main__':

    print " ".join(sys.argv)

    # Command-line input
    parser = argparse.ArgumentParser(description="Generate a guessed force field")
    parser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    parser.add_argument('-i', '--itp', help="an ITP with the all-atom structure")
    parser.add_argument('-nb', '--nonbonded', help="an ITP with the non-bonded parameters")
    args = parser.parse_args()

    # Load the force field definitions
    ff = elbalib.Elba()
    ff.load(args.xml)

    # Read the all atom file
    aaitp = gmxlib.SimpleItp()
    aaitp.read_itp(args.itp)
    aaitp.read_atomtypes(args.nonbonded)

    # We will make a copy of the force field molecule that we can modify
    ffmol = ff.find_molecule(aaitp.name)
    ffmol_new = copy.deepcopy(ffmol)

    # CHARGES
    print "[ Charges ]"
    print "%4s %5s %5s"%("name","new","old")
    for bn in ffmol.beadnames:
        b1 = ffmol_new.beads[bn]
        if b1.beadtype.charge not in [0.0,"x"] :
            atoms = [aaitp.find_atom(a) for a in b1.chargegroup]
            b1.beadtype.charge = 0
            for a in atoms :
                b1.beadtype.charge += a.charge
                b1.beadtype.charge += sum([a2.charge if a2.type.atomnum == 1 else 0.0 for a2 in aaitp.bonded_to(a)])

            b2 = ffmol.beads[bn]
            print "%4s %5.2f %s"%(bn,b1.beadtype.charge,"%5.2f"%b2.beadtype.charge if b2.beadtype.charge != "x" else "%5s"%b2.beadtype.charge)


    # LJ parameters
    c12_elbaw,c6_elbaw=elbalib.ljconvert(epsilon=0.550,  sigma=3.050)
    c12_fac = c12_elbaw/600000.0
    c6_fac  = c6_elbaw/610.0
    print "\n[ LJ parameters ]",c12_fac,c6_fac
    for bn in ffmol.beadnames:
        b1 = ffmol_new.beads[bn]
        c12_old,c6_old = elbalib.ljconvert(sigma=b2.beadtype.sigma,epsilon=b2.beadtype.epsilon)
        #

        atoms = []
        for a in b1.mapping:
            atoms.append(aaitp.find_atom(a))
            for a2 in aaitp.bonded_to(atoms[-1]):
                if a2.type.atomnum == 1 : atoms.append(a2)
        c6,c12 = 0.0,0.0
        for ai,a in enumerate(atoms) :
            #print a.name,a.type.sigma,
            ac12,ac6 = elbalib.ljconvert(sigma=a.type.sigma,epsilon=a.type.epsilon)
            c6 += ac6
            c12 += ac12
        #c12 *= c12_fac
        #c6 *= c6_fac
        b1.beadtype.epsilon,b1.beadtype.sigma = elbalib.ljconvert(c6=c6,c12=c12)

        print "%4s %5.3f %5.3f %10.3e %10.3e"%(bn,b1.beadtype.epsilon,b1.beadtype.sigma,c12,c6),
        b2 = ffmol.beads[bn]
        print " %5.3f %5.3f %10.3e %10.3e"%(b2.beadtype.epsilon,b2.beadtype.sigma,c12_old,c6_old)

    print
