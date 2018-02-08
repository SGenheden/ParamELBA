# Author: Samuel Genheden samuel.genheden@gmail.com

import argparse
import os
import sys

import MDAnalysis as md
import numpy as np

import elbalib

def generate_lmpinput(structure,ff,postfix) :

    # Cach up ff molecules, convert coordinates and count atoms and topology elements
    xyz = []
    ffmols = []
    natoms = 0
    nbonds = 0
    nangles = 0
    nwat = 0
    print len(structure.selectAtoms("all").residues)
    for residue in structure.selectAtoms("all").residues :
        if residue.name not in ff.molecules :
            raise Exception("Could not find %s in force field"%residue.name)
        if residue.name == "SOL" : nwat += 1
        ffmols.append(ff.molecules[residue.name])
        trmat = ffmols[-1].transmat([a.name.lower() for a in residue])
        xyz.append(np.dot(trmat,residue.get_positions()))
        natoms += len(ffmols[-1].beads)
        nbonds += len(ffmols[-1].bonds)
        nangles += len(ffmols[-1].angles) + len(ffmols[-1].dipolerestraints)

    print nwat

    with open("data."+postfix,"w") as f :

        f.write("Created by generate_lmpinput.py from %s\n"%postfix)
        f.write("\n")
        f.write("%d atoms\n"%natoms)
        f.write("%d bonds\n"%nbonds)
        f.write("%d angles\n"%nangles)
        f.write("\n")
        f.write("%d atom types\n"%len(ff.beadtypes))
        f.write("%d bond types\n"%len(ff.bondtypes))
        f.write("%d angle types\n"%(len(ff.angletypes)+len(ff.dipoletypes)))
        f.write("\n")

        box = structure.dimensions[:3]
        f.write("%15.5f %15.5f xlo xhi\n"%(0.0,box[0]))
        f.write("%15.5f %15.5f ylo yhi\n"%(0.0,box[1]))
        f.write("%15.5f %15.5f zlo zhi\n"%(0.0,box[2]))
        f.write("\n")

        f.write("Atoms\n\n")
        natoms = 0
        for i,(mol,mpos) in enumerate(zip(ffmols,xyz),1) :
            for bn,apos in zip(mol.beadnames,mpos) :
                natoms += 1
                f.write("%s\n"%mol.beads[bn].generate_input(natoms,i,apos))

        if nbonds > 0 :
            f.write("\nBonds\n\n")
            natoms = 1
            nbonds = 0
            for mol in ffmols :
                for bond in mol.bonds :
                    nbonds += 1
                    f.write("%6d %2d %s\n"%(nbonds,bond.contype.id," ".join(["%6d"%(i+natoms) for i in bond.indices(mol)])))
                natoms += len(mol.beads)

        if nangles > 0 :
            f.write("\nAngles\n\n")
            natoms = 1
            nangles = 0
            for mol in ffmols :
                for angle in mol.angles :
                    nangles += 1
                    f.write("%6d %2d %s\n"%(nangles,angle.contype.id," ".join(["%6d"%(i+natoms) for i in angle.indices(mol)])))
                for diprest in mol.dipolerestraints :
                    nangles += 1
                    idxs = diprest.indices(mol)+natoms
                    f.write("%6d %2d %s %6d\n"%(nangles,diprest.contype.id," ".join(["%6d"%i for i in idxs]),idxs[-1]))
                natoms += len(mol.beads)

        f.write("\n\n")

    with open("bonded."+postfix,"w") as f :

        f.write("# harmonic bond coefficients:\n")
        f.write("#         bondType   K     r0\n")
        btns = sorted(ff.bondtypes)
        for b in btns:
            typ = ff.bondtypes[b]
            f.write("bond_coeff  %5d %8.3f %8.3f\n"%(b,typ.k,typ.r0))

        f.write("\n# angle coefficients:\n")
        f.write("#         angleType                  K    theta0\n")
        atns = sorted(ff.angletypes)
        for a in atns:
            typ = ff.angletypes[a]
            f.write("angle_coeff  %5d cosine/squared %8.3f %8.3f\n"%(a,typ.k,typ.theta0))
        atns = sorted(ff.dipoletypes)
        for a in atns:
            typ = ff.dipoletypes[a]
            f.write("angle_coeff  %5d dipole %8.3f %8.3f\n"%(a,typ.k,typ.gamma0))

def modify_bonded(ff,resname,datafile,modifiers=None) :

    allmod = {tag:True for tag in "bond_k bond_equil ang_k ang_equil".split()}
    if modifiers is None:
        modifiers = {k:v for k,v in allmod.iteritems()}
    else:
        for k in allmod:
            if k not in modifiers:
                modifiers[k] = True
    print modifiers
    ffmol = ff.molecules[resname]

    data = []
    with open(datafile,"r") as f :
        for line in f.readlines() :
            data.append(line.strip().split())
    data = np.array(data,dtype=float)

    di = 0
    ff.bondtypes = {}
    for bond in ffmol.bonds :
        k = data[di,0] if modifiers["bond_k"] else bond.contype.k
        r0 =  data[di,1] if modifiers["bond_equil"] else bond.contype.r0
        bond.contype = elbalib.BondType()
        bond.contype.set(id=di+1,kf=1,k=k,r0=r0)
        ff.bondtypes[di+1]=bond.contype
        di += 1

    ai = 0
    prev = len(ff.angletypes)
    ff.angletypes = {}
    for angle in ffmol.angles :
        k =  data[di,0] if modifiers["ang_k"] else angle.contype.k
        theta0 = data[di,1] if modifiers["ang_equil"] else angle.contype.theta0
        angle.contype = elbalib.AngleType()
        angle.contype.set(id=ai+1,kf=1,k=k,theta0=theta0)
        ff.angletypes[ai+1]=angle.contype
        di += 1
        ai += 1

    offset = len(ff.angletypes)-prev
    lst = {}
    for dt in ff.dipoletypes :
        ff.dipoletypes[dt].id += offset
        lst[ff.dipoletypes[dt].id] = ff.dipoletypes[dt]
    ff.dipoletypes = lst


if __name__ == '__main__':

    print " ".join(sys.argv)

    argparser = argparse.ArgumentParser(description="Program the generate LAMMPS input files")
    argparser.add_argument('-s', '--struct', help="a structure file")
    argparser.add_argument('-x', '--xml', help="an XML file with force field definitions")
    argparser.add_argument('--bonded',nargs=2,help="replace bonded parameters")
    argparser.add_argument('--modifiers',nargs="+",help="modifiers for replacement")
    args = argparser.parse_args()

    struct = md.Universe(args.struct)
    ff = elbalib.Elba()
    ff.load(args.xml)

    if args.bonded :
        modifiers = None
        if args.modifiers is not None:
            modifiers={k:False for k in args.modifiers}
        modify_bonded(ff,args.bonded[0],args.bonded[1],modifiers)

    postfix = os.path.splitext(os.path.basename(args.struct))[0]
    generate_lmpinput(struct,ff,postfix)
