# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse
import os
import sys

import numpy as np

import elbalib

def generate_gmxinput(ffmol) :

    print "[bonds]"
    print "; i j   funct   length  force.c."
    for bond in ffmol.bonds :
        print "%2d %2d"%tuple(bond.indices(ffmol)+1),
        print "  1 %8.2f %8.2f"%(bond.contype.r0*0.1,bond.contype.k*4.184*2*100)

    print "[angles]"
    print "; i j k         funct   angle   force.c."
    for angle in ffmol.angles :
        print "%2d %2d %2d"%tuple(angle.indices(ffmol)+1),
        print "   2 %8.2f %8.2f"%(angle.contype.theta0,angle.contype.k*4.184*2)

def modify_bonded(ff,resname,datafile) :

    ffmol = ff.molecules[resname]

    data = []
    with open(datafile,"r") as f :
        for line in f.readlines() :
            data.append(line.strip().split())
    data = np.array(data,dtype=float)

    di = 0
    ff.bondtypes = {}
    for bond in ffmol.bonds :
        bond.contype = elbalib.BondType()
        bond.contype.set(id=di+1,kf=2,k=data[di,0],r0=data[di,1])
        ff.bondtypes[di+1]=bond.contype
        di += 1

    ai = 0
    prev = len(ff.angletypes)
    ff.angletypes = {}
    for angle in ffmol.angles :
        angle.contype = elbalib.AngleType()
        angle.contype.set(id=ai+1,kf=2,k=data[di,0],theta0=data[di,1])
        ff.angletypes[ai+1]=angle.contype
        di += 1
        ai += 1


if __name__ == '__main__':

    print " ".join(sys.argv)

    argparser = argparse.ArgumentParser(description="Program the generate LAMMPS input files")
    argparser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    argparser.add_argument('-b','--bonded',nargs=2,help="replace bonded parameters")
    args = argparser.parse_args()

    ff = elbalib.Elba()
    ff.load(args.xml)

    modify_bonded(ff,args.bonded[0],args.bonded[1])
    generate_gmxinput(ff.molecules[args.bonded[0]])
