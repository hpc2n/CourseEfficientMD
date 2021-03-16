#!/usr/bin/env python
"""
charmmgui2dms.py: convert a charmm system to a Desmond dms file or a LAMMPS data file

Correspondance: qiyifei@gmail.com or wonpil@lehigh.edu

Last update: May 12, 2017

"""


import sys
import os
import sqlite3
import glob
import argparse
from argparse import RawTextHelpFormatter
import numpy as np
import math
import itertools
import datetime
import getpass


class atom:
    atom_mass = """1 H 1.0079
    2 He 4.0026
    3 Li 6.941
    4 Be 9.0122
    5 B 10.811
    6 C 12.0107
    7 N 14.0067
    8 O 15.9994
    9 F 18.9984
    10 Ne 20.1797
    11 Na 22.9897
    12 Mg 24.305
    13 Al 26.9815
    14 Si 28.0855
    15 P 30.9738
    16 S 32.065
    17 Cl 35.453
    19 K 39.0983
    18 Ar 39.948
    20 Ca 40.078
    21 Sc 44.9559
    22 Ti 47.867
    23 V 50.9415
    24 Cr 51.9961
    25 Mn 54.938
    26 Fe 55.845
    28 Ni 58.6934
    27 Co 58.9332
    29 Cu 63.546
    30 Zn 65.39
    31 Ga 69.723
    32 Ge 72.64
    33 As 74.9216
    34 Se 78.96
    35 Br 79.904
    36 Kr 83.8
    37 Rb 85.4678
    38 Sr 87.62
    39 Y 88.9059
    40 Zr 91.224
    41 Nb 92.9064
    42 Mo 95.94
    43 Tc 98
    44 Ru 101.07
    45 Rh 102.9055
    46 Pd 106.42
    47 Ag 107.8682
    48 Cd 112.411
    49 In 114.818
    50 Sn 118.71
    51 Sb 121.76
    53 I 126.9045
    52 Te 127.6
    54 Xe 131.293
    55 Cs 132.9055
    56 Ba 137.327
    57 La 138.9055
    58 Ce 140.116
    59 Pr 140.9077
    60 Nd 144.24
    61 Pm 145
    62 Sm 150.36
    63 Eu 151.964
    64 Gd 157.25
    65 Tb 158.9253
    66 Dy 162.5
    67 Ho 164.9303
    68 Er 167.259
    69 Tm 168.9342
    70 Yb 173.04
    71 Lu 174.967
    72 Hf 178.49
    73 Ta 180.9479
    74 W 183.84
    75 Re 186.207
    76 Os 190.23
    77 Ir 192.217
    78 Pt 195.078
    79 Au 196.9665
    80 Hg 200.59
    81 Tl 204.3833
    82 Pb 207.2
    83 Bi 208.9804
    84 Po 209
    85 At 210
    86 Rn 222
    87 Fr 223
    88 Ra 226
    89 Ac 227
    91 Pa 231.0359
    90 Th 232.0381
    93 Np 237
    92 U 238.0289
    95 Am 243
    94 Pu 244
    96 Cm 247
    97 Bk 247
    98 Cf 251
    99 Es 252
    100 Fm 257
    101 Md 258
    102 No 259
    104 Rf 261
    103 Lr 262
    105 Db 262
    107 Bh 264
    106 Sg 266
    109 Mt 268
    111 Rg 272
    108 Hs 277"""

    temp = atom_mass.split()
    anum = []
    amass = []
    for i in range(0, len(temp), 3):
        anum.append(int(temp[i]))
        amass.append(float(temp[i+2]))
    amass = np.array(amass, dtype=np.float)

    def __init__(self, atomid, segid, resid, resname, atomname, atomtype, charge, mass):
        self.atomid = int(atomid)
        self.segid  = segid
        self.resid  = resid
        self.resname = resname
        self.atomname = atomname
        #self.atomname = atomname.replace("'","")
        self.atomtype = atomtype
        self.charge = float(charge)
        self.mass = float(mass)
        self.coor = [9999,9999,9999]
        self.atomstr = "-".join([str(self.atomid), self.segid, str(self.resid),\
                self.resname, self.atomname, self.atomtype, str(self.charge)])
        self.atomstr = "-".join([self.segid, str(self.resid),\
                self.resname, self.atomname, self.atomtype, str(self.charge)])
        self.anum = self.get_anum(self.mass, self.atomtype)
        self.nonbond_parid = -1

    def get_info(self):
        return " ".join([self.segid, self.resname, self.resid, self.atomname])
    def set_coor(self, x, y, z):
        self.coor = [x, y, z]
    def get_coor(self):
        return self.coor

    def get_anum(self, atommass, atomtype):
        idx = (np.abs(atom.amass-atommass)).argmin()
        return atom.anum[idx]

class psf:
    def __init__(self, psffile):
        self.atom = []
        self.bond = []
        self.angle = []
        self.dihe = []
        self.impr = []
        self.cmap = []
        self.read_psf(psffile)

    def stat(self):
        print "read from psf file %s"%self.psffile
        print len(self.atom),"atoms"
        print len(self.bond),"bonds(not include TIP3-H -- TIP3-H bond in psf)"
        print len(self.angle),"angles"
        print len(self.dihe),"dihes"
        print len(self.impr),"impropers"
        print len(self.cmap),"CMAPs"

    def read_psf(self,psffile):
        fp = open(psffile, "r")
        for line in fp.readlines():
            line = line.strip("\n")
            if len(line) == 0:
                continue
            if line[0] == "*":
                continue
            linee = line.split()
            if len(linee) == 0:
                continue
            if linee[0] == "PSF":
                continue
            n = len(linee)
            if n > 1 and linee[1] == "!NTITLE":
                tag = "!NTITLE"
                continue
            if n > 1 and linee[1] in ["!NATOM", "!NBOND:", "!NTHETA:", \
                    "!NPHI:","!NIMPHI:", "!MOLNT", "!NCRTERM:", "!NNB", "!NDON:", "!NACC:"]:
                tag = linee[1]
                continue
            if len(linee) > 2 and linee[2] in ["!NGRP", "!NUMLP"]:
                tag = linee[2]
                continue
            if tag == "!NATOM":
            #  1 PROA     417      LYS      N        NH3     -0.300000       14.0070           0   0.00000     -0.301140E-02
                if linee[5].isdigit():
                    raise Exception("atom type must not be a number: %s"%line)
                atm = atom(linee[0], linee[1], linee[2], linee[3], linee[4], linee[5], linee[6], linee[7])
                self.atom.append(atm)
                continue
            if tag == "!NTHETA:":
                for i in range(0, len(linee), 3):
                    id1 = int(linee[i])
                    id2 = int(linee[i+1])
                    id3 = int(linee[i+2])
                    self.angle.append((id1, id2, id3))
            if tag == "!NPHI:":
                for i in range(0, len(linee), 4):
                    id1 = int(linee[i])
                    id2 = int(linee[i+1])
                    id3 = int(linee[i+2])
                    id4 = int(linee[i+3])
                    self.dihe.append((id1, id2, id3, id4))

            if tag == "!NIMPHI:":
                for i in range(0, len(linee), 4):
                    id1 = int(linee[i])
                    id2 = int(linee[i+1])
                    id3 = int(linee[i+2])
                    id4 = int(linee[i+3])
                    self.impr.append((id1, id2, id3, id4))

            if tag == "!NBOND:":
                for i in range(0, len(linee), 2):
                    id1 = int(linee[i])
                    id2 = int(linee[i+1])
                    if self.atom[id1-1].resname == "TIP3" and self.atom[id1-1].atomname[0]=="H"\
                            and self.atom[id2-1].atomname[0]=="H":
                                #bond between TIP3 H-H removed
                                continue
                    else:
                    #self.bond.add((id1, id2))
                        self.bond.append((id1, id2))
                        
            if tag == "!NCRTERM:":
                for i in range(0, len(linee), 8):
                    self.cmap.append( (int(linee[i]), int(linee[i+1]), int(linee[i+2]), int(linee[i+3]),\
                            int(linee[i+4]), int(linee[i+5]), int(linee[i+6]), int(linee[i+7])))
        self.atomnum = len(self.atom)
        self.bondnum = len(self.bond)
        self.anglenum = len(self.angle)
        self.dihenum = len(self.dihe)
        self.imprnum = len(self.impr)
        self.cmapnum = len(self.cmap)
        self.psffile = psffile
        print "read %d atoms %d bonds %d angle %d dihe %d impr  %d cmap from %s"%(self.atomnum, self.bondnum, self.anglenum, self.dihenum, self.imprnum, self.cmapnum, psffile)

    def get_uniq_dihedral_pairs(self):
        dih_pairs = []
        for i in self.dihe:
            types = [ self.atom[j-1].atomtype for j in i ]
            type1 = " ".join(types)
            if not type1 in dih_pairs:
                dih_pairs.append(type1)
        return dih_pairs

    def get_uniq_impr_pairs(self):
        impr_pairs = []
        for i in self.impr:
            types = [ self.atom[j-1].atomtype for j in i ]
            type1 = " ".join(types)
            if not type1 in impr_pairs:
                impr_pairs.append(type1)
        return impr_pairs

    def get_parameter(self, prm):
        #nonbond type
        self.nonbond_parid = []
        nonbond_parid_dict = {}
        for i in range(0, self.atomnum):
            atomtype = self.atom[i].atomtype
            if atomtype in nonbond_parid_dict:
                id = nonbond_parid_dict[atomtype]
            else:
                id =prm.get_nonbond_par_index(atomtype)
                nonbond_parid_dict[atomtype] = id
            self.nonbond_parid.append(id)
        
        #bond type
        self.bond_parid = []
        bond_parid_dict = {}
        for i in range(0, self.bondnum):
            type1 = self.atom[self.bond[i][0]-1].atomtype
            type2 = self.atom[self.bond[i][1]-1].atomtype
            bondtype = type1 + " " + type2
            if bondtype in bond_parid_dict:
                id = bond_parid_dict[bondtype]
            else:
                id = prm.get_bond_par_index(type1, type2)
                bond_parid_dict[bondtype] = id
                bondtype = type2 + " " + type1
                bond_parid_dict[bondtype] = id
            self.bond_parid.append(id)
       
       #angle type
        self.angle_parid = []
        angle_parid_dict = {}
        for i in range(0, self.anglenum):
            type1 = self.atom[self.angle[i][0]-1].atomtype
            type2 = self.atom[self.angle[i][1]-1].atomtype
            type3 = self.atom[self.angle[i][2]-1].atomtype
            angletype = type1 + " " + type2 + " " + type3
            if angletype in angle_parid_dict:
                id = angle_parid_dict[angletype]
            else:
                id = prm.get_angle_par_index(type1, type2, type3)
                angle_parid_dict[angletype] = id
                angletype = type3 + " " + type2 + " " + type1
                angle_parid_dict[angletype] = id
            self.angle_parid.append(id)

       #dihedral type
        self.dihe_parid = []
        dihe_parid_dict = {}
        for i in range(0, self.dihenum):
            type1 = self.atom[self.dihe[i][0]-1].atomtype
            type2 = self.atom[self.dihe[i][1]-1].atomtype
            type3 = self.atom[self.dihe[i][2]-1].atomtype
            type4 = self.atom[self.dihe[i][3]-1].atomtype
            dihetype = " ".join([type1, type2, type3, type4])
            if dihetype in dihe_parid_dict:
                id = dihe_parid_dict[dihetype]
            else:
                id = prm.get_dihe_par_index(type1, type2, type3, type4)
                dihe_parid_dict[dihetype] = id
                dihetype = " ".join([type4, type3, type2, type1])
                dihe_parid_dict[dihetype] = id
            self.dihe_parid.append(id)

        #impr type
        self.impr_parid = []
        impr_parid_dict= {}
        for i in range(0, self.imprnum):
            type1 = self.atom[self.impr[i][0]-1].atomtype
            type2 = self.atom[self.impr[i][1]-1].atomtype
            type3 = self.atom[self.impr[i][2]-1].atomtype
            type4 = self.atom[self.impr[i][3]-1].atomtype
            imprtype = " ".join([type1, type2, type3, type4])
            if imprtype in impr_parid_dict:
                id = impr_parid_dict[imprtype]
            else:
                id = prm.get_impr_par_index(type1, type2, type3, type4)
                impr_parid_dict[imprtype] = id
                imprtype = " ".join([type4, type3, type2, type1])
                impr_parid_dict[imprtype] = id
            self.impr_parid.append(id)

        #CMAP type
        self.cmap_parid = []
        for i in range(0, self.cmapnum):
            type = [self.atom[atomid-1].atomtype for atomid in self.cmap[i]]
            cmaptype = " ".join(type)
            id = prm.get_cmap_par_index(type)
            self.cmap_parid.append(id)

class crd:
    def __init__(self, crdfile):
        self.atom = []
        self.read_crd(crdfile)

    def read_crd(self, crdfile):
        fp = open(crdfile, "r")
        natom = 0
        for line in fp.readlines():
            line = line.strip("\n")
            if line[0] == "*":
                continue
            linee= line.split()
            if len(linee) == 2:
                natom = int(linee[0])
                continue
            if len(linee) == 1:
                natom = int(linee[0])
                continue
            #         1         1  SER       N               6.5056730000       50.0926600000        5.8024520000  MSP1      55              0.0000000000
            atomid, residindex, resname, atomname, x, y, z, segid, resid, temp = linee 
            atm = atom(atomid, segid, resid, resname, atomname, atomtype="", charge=999, mass=-1)
            atm.set_coor(float(x), float(y), float(z))
            self.atom.append(atm)

        print "read in",len(self.atom), "atoms from", crdfile
        if len(self.atom) != natom:
            print "number of atom does not match crd head information"
        self.atomnum = len(self.atom)

    def add_posk(self, posk):
        self.posk = posk[:]

    def get_coor(self, i):
        return self.atom[i].get_coor()


class prm_atom:
    def __init__(self, name, mass, memo = ""):
        self.name = name
        self.mass = mass
        self.memo = memo.replace("'"," ")
    def __eq__(self, other):
        if self.name == other.name:
            return True
        else:
            return False

class prm_bond:
    def __init__(self, atoms, kb, b0, memo = ""):
        self.atom = atoms[:]
        self.kb = kb
        self.b0 = b0
        self.memo = memo.replace("'", " ")
    def __eq__(self, other):
        if self.name == other.name or self.name == other.name[::-1]:
            return True
        else:
            return False

class prm_angle:
    def __init__(self, atoms, Ktheta, Theta0, Kub=-1, S0=-1, memo=""):
        self.atom   = atoms[:]
        self.Ktheta = Ktheta
        self.Theta0 = Theta0
        self.Kub    = Kub
        self.S0     = S0
        self.memo = memo.replace("'", " ")
    def __eq__(self, other):
        if self.atom == other.atom or self.atom == other.atom[::-1]:
            return True
        else:
            return False

class prm_dihe:
    """V_dihe_anton = c0 + sum(cn*cos(n*phi-delta))
    V_dihe_charmm = sum(kn * (1 + cos(n*phi - delta)))
    so:
    c0 = sum(kn)
    if delta == 180 in charmm, convert it to 0.0, and turn k to -k in self.fc[1,2,3,4,5,6]
    """

    def __init__(self, names, kchi, n, delta, memo = ""):
        self.names = names[:]
        self.dihe_orig =[ [n, kchi, delta] ]
        #one dihedral may have different deltas
        #put them in to different dihedral terms
        #delta = 180 and 0 belong to the same dihedral term
        if delta == 180:
            delta_str = str(0)
        else:
            delta_str = str(int(delta))
        self.namestr = "_".join(self.names + [delta_str])
        self.added_mult = [n]
        self.memo = memo.replace("'", " ")
        self.lammps_weight = -1
        self.generate_desmond()


    def generate_lammps_coef(self, weight):
        #      0.2          3          0          1  # HA   CT3  CT1  OH1 
        prmline = []
        names = " ".join(self.names)
        for i in range(0, len(self.dihe_orig)):
            n = self.dihe_orig[i]
            #only the first n has weight, others have 0 weight
            if i == 0:
                w = weight
            else:
                w = 0
            #k n delta weight
            line = "%12.4f %8d %12d %3.2f #%s"%(n[1], n[0], n[2], w, names )
            prmline.append(line)
        return prmline

    def generate_desmond(self):
        self.fc = [0.0] * 7
        self.delta = "none"
        for i in self.dihe_orig:
            #[n, kchi, delta]
            self.fc[0] += i[1]

        for i in self.dihe_orig:
            mult = i[0]
            k = i[1]
            d = i[2]
            if d == 180:
                delta = 0
                k = -k
            else:
                delta = d
            if self.delta == "none":
                self.delta = delta
            else:
                if self.delta != delta:
                    raise Exception("generate desmond dihe prm error")
            self.fc[mult] = k

    def add_mult(self, names, kchi, n, delta):
        if delta == 180:
            delta_str = str(0)
        else:
            delta_str = str(int(delta))
        namestr = "_".join(names + [delta_str])
        if namestr != self.namestr:
            raise Exception("name does not match trying add dihe %s to %s"%(namestr, self.namestr))
 
        if [n, kchi, delta] in self.dihe_orig:
            print "warnning duplicate dihe term %s, will skip"%(" ".join(self.names)), self.dihe_orig
            print "try to add", names, n, kchi, delta
            return

        tag = False
        for i in self.dihe_orig:
            # [n, kchi, delta]
            if n == i[0]:
                tag = True
                k = i[1]
                d = i[2]
                index = self.dihe_orig.index(i)
                self.dihe_orig[index][1] = kchi
                self.dihe_orig[index][2] = delta
                print "warnning overwriting dihe term %s, multi %d, chi %f delta %f -> chi %f delta %f"%(self.namestr, n,\
                        k, d, kchi, delta)
        if not tag:
            self.dihe_orig.append([n, kchi, delta])
            self.dihe_orig.sort(key=lambda x: x[0])
            self.added_mult.append(n)
        self.generate_desmond()


class prm_impr:
    def __init__(self, atoms, kpsi, psi0, memo = ""):
        self.atoms = atoms[:]
        self.kpsi = kpsi
        self.psi0 = psi0
        self.memo = memo.replace("'", " ")
    def __eq__(self, other):
        if self.atoms == other.atoms[::-1]:
            return True
        else:
            return False

class prm_nonbond:
    def __init__(self, name, epsilon, sigma, epsilon_14, sigma_14, memo = ""):
        self.name = name
        self.epsilon = epsilon
        if self.epsilon < 0:
            raise Exception("epsilon in nonbond should > 0 in Anton format")
        self.sigma = sigma
        self.epsilon_14 = epsilon_14
        #sigma_14 < 0 means no 1_4 nonbond parameters
        self.sigma_14 = sigma_14
        self.memo = memo.replace("'", " ")

    def __eq__(self, other):
        if self.name == other.name:
            return True
        else:
            return False

class prm_cmap:
    def __init__(self, id, atomlist, ntorsion):
        if not len(atomlist) == 8:
            raise Exception("CMAP must have 8 atoms, the input is %s"%("-".join(namelist)))
        if not ntorsion == 24:
            raise Exception("only CMAP 24 is supported, the input is %d"%ntorsion)
        self.id = id
        self.atom = atomlist[:]
        self.ntorsion = ntorsion
        self.psi=np.linspace(-180, 180, ntorsion, endpoint=False)
        self.cmap = []
        self.phi=np.linspace(-180, 180, ntorsion, endpoint=False)

    def add_torsion(self, values):
        #self.phi.append(phi)
        #if not len(values) == self.ntorsion:
        #    raise Exception("number of cmap %d values does not match ntorsion %d"%(len(values),self.ntorsion))
        self.cmap.extend(values[:])
    #! phi = -180.0
    #     0.126790      0.768700      0.971260      1.250970      2.121010
    #     2.720430      2.089440      1.789790      0.780870     -0.688474
    #     1.001130     -2.200520     -4.827670     -4.821447     -4.913223
    #    -3.591106     -2.766446     -2.784200     -2.454589     -2.346991
    #    -2.335350     -1.522656     -0.951542     -0.036650 

    def __eq__(self, other):
        if self.atom == other.atom:
            return True
        else:
            return False

    def write(self, fp):

        fp.write("#CMAP for %s; id=%d\n"%(" ".join(self.atom), self.id))
        for i in range (0, len(self.phi)):
            fp.write("#phi = %f\n"%self.phi[i])
            count = 1
            for j in range(i*24, i*24+24):
                fp.write(" %12.6f"%(self.cmap[j]))
                if count % 5 == 0:
                    fp.write("\n")
                count += 1
            fp.write("\n\n")
        fp.write("\n")


class prm_nbfix:
    def __init__(self, type1, type2, emin, rmin, emin14, rmin14, memo=""):
        self.type = [type1, type2]
        self.emin = float(emin)
        self.rmin = float(rmin)
        self.emin14 = float(emin14)
        self.rmin14 = float(rmin14)
        self.memo = memo.replace("'", " ")
        #E_charmm = emin(rmin^12/r^12 - 2 * rmin^6/r^6)
        #E_desmond = 4epsilon*(sigma^12/r^12 - sigma^6/r^6)
        self.sigma = (self.rmin**6/2)**(1.0/6.0)
        if self.rmin14 > 0:
            self.sigma14 = (self.rmin14**6/2)**(1.0/6.0)
        else:
            self.sigma14 = -1
        self.epsilon_desmond = -self.emin
        self.epsilon14_desmond = -self.emin14

    def __eq__(self, other):
        if self.type == other.type:
            return True
        else:
            return False

def isfloat(str):
    try:
        temp = float(str)
        return True
    except:
        return False

class prm:
    def __init__(self, prmfile = ""):
        self.atom    = []
        self.atom_dict = {}
        self.bond    = []
        self.bond_dict = {}
        self.angle   = []
        self.angle_dict = {}
        self.dihe    = []
        self.dihe_dict = {}
        self.dihe_dict["wildcard"] = {}
        self.impr    = []
        self.impr_dict = {}
        self.impr_dict["wildcard"] = {}
        self.cmap    = []
        self.cmap_dict = {}
        self.nonbond = []
        self.nonbond_dict = {} 
        self.nbfix = []
        self.nbfix_dict = {}
        self.prmfiles = []

        if prmfile != "":
            self.read_prm(prmfile)
    
    def get_mass_par_index(self, type):
        type = type.upper()
        if not type in self.atom_dict:
            raise Exception("mass parameter not found for %s"%type)
        return self.atom_dict[type]

    def get_mass(self, type):
        atomid = self.get_mass_par_index(type)
        return self.atom[atomid].mass

    def get_nonbond_par_index(self, type):
        type = type.upper()
        if not type in self.nonbond_dict:
            raise Exception("nonbond parameter not found for %s"%type)
        return self.nonbond_dict[type]

    def get_bond_par_index(self, type1, type2):
        type = type1.upper()+" "+type2.upper()
        if type1.upper() in ["X","*"] or type2.upper() == ["X", "*"]:
            raise Exception("wildcard not support in prm.get_bond_par_index for type %s"%type)
        if not type in self.bond_dict:
            raise Exception("bond parameter not found for %s"%type)
        #this is the index of the bond para in self.bond
        return self.bond_dict[type]

    def get_angle_par_index(self, type1, type2, type3):
        type1 = type1.upper()
        type2 = type2.upper()
        type3 = type3.upper()
        types = [type1, type2, type3]
        type = " ".join(types)
        if "X" in types or "*" in types:
            raise Exception("wildcard not support in prm.get_angle_par_index for type %s"%type)
        if not type in self.angle_dict:
            raise Exception("angle parameter not found for %s"%type)
        #this is the index of the angle para in self.angle
        return self.angle_dict[type]

    def wildcard_match(self, name1, name2):
        match = 0
        name1 = name1.replace("*", "X").upper().split()
        name2 = name2.replace("*", "X").upper().split()
        assert len(name1) == len(name2)
        for k in range(0, len(name1)):
            if name1[k] == "X" or name2[k] == "X":
                match += 1
                continue
            elif name1[k] == name2[k]:
                match += 1

        if match == len(name1):
            return True
        else:
            return False

    def get_dihe_par_index(self, type1, type2, type3, type4):
        type = [type1.upper(), type2.upper(), type3.upper(), type4.upper()]
        type1 = " ".join(type)
        if type1 in self.dihe_dict:
            return self.dihe_dict[type1]
        else:
            for wildcard_dihe in self.dihe_dict["wildcard"]:
                if self.wildcard_match(type1, wildcard_dihe):
                    return self.dihe_dict["wildcard"][wildcard_dihe]
        for i in self.dihe_dict["wildcard"]:
            print i, self.dihe_dict["wildcard"][i]
        raise Exception("dihedarl param for %s not found"%type1)

    def get_impr_par_index(self, type1, type2, type3, type4):
        #allowing reverse order match
        type = [type1.upper(), type2.upper(), type3.upper(), type4.upper()]
        type1 = " ".join(type)
        type2 = " ".join(type[::-1])
        if type1 in self.impr_dict:
            return self.impr_dict[type1]
        elif type2 in self.impr_dict:
            return self.impr_dict[type2]
        else:
            for wildcard_impr in self.impr_dict["wildcard"]:
                if self.wildcard_match(type1, wildcard_impr):
                    return self.impr_dict["wildcard"][wildcard_impr]
                if self.wildcard_match(type2, wildcard_impr):
                    return self.impr_dict["wildcard"][wildcard_impr]

        for i in self.impr_dict:
            if i != "wildcard": 
                print i, self.impr_dict[i]
        for i in self.impr_dict["wildcard"]:
            print i, self.impr_dict["wildcard"][i]

        raise Exception("impr param for %s not found"%type1)

    def get_cmap_par_index(self, atomtypes):
        assert len(atomtypes) ==  8
        type = " ".join(atomtypes).upper()
        if not type in self.cmap_dict:
            raise Exception("CMAP parameter not found for %s"%type)
        return self.cmap_dict[type]

    def get_nbfix_par_index(self, type1, type2):
        type = type1 + " " + type2
        if type.upper() in self.nbfix_dict:
            return self.nbfix_dict[type]
        else:
            return -1

    def stat(self):
        print "charmm FF files are:",self.prmfiles
        print "read in"
        print "%d atom types"%len(self.atom)
        print "%d bond par"%len(self.bond)
        print "%d angle par"%len(self.angle)
        print "%d dihe par"%len(self.dihe)
        print "%d impr par"%len(self.impr)
        print "%d cmap par"%len(self.cmap)
        print "%d nonbond par"%len(self.nonbond)
        print "%d nbfix par"%len(self.nbfix)
        print

    def read_prm(self, prmfile):
        fp = open (prmfile, "r")
        self.prmfiles.append(prmfile)
        tags = ["ATOMS", "BONDS", "ANGLES", "DIHEDRALS", "IMPROPER", "NONBONDED", "CMAP", "NBFIX", "IMPROPERS"]
        curtag = ""
        for line in fp.readlines():
            line = line.strip("\n")
            line = line.lstrip(" ")
            if len(line) == 0 or line[0] == "!" or line[0] == "*" :
                continue
            linee = line.split()
            if len(linee) == 0:
                continue
            if linee[0].upper() in ["READ","END", "CUTNB", "HBOND", "RETURN","SET", \
                    "IF","DONOR","ACCEPTOR","IC","RESI","PRES","PATCH","GROUP","PATCHING"\
                    ,"ATOM","BOND","DOUBLE","IMPR","DELETE","DECL","DEFAULT","TRIPLE"]:
                continue

            if len(linee[0]) >= 4 and linee[0].upper()[0:4] in ["CUTN", "HBON","RETU","DONO","ACCE"\
                    ,"PATC","GROU","PATC","DOUB","DELE","DECL","DEFA","BILD","TRIP"]:
                continue

            if linee[0].upper() in tags:
                if linee[0].upper() == "CMAP":
                    #print linee
                    try:
                        id = linee.index("!")
                    except:
                        id = -1
                    if len(linee[0:id]) > 1:
                        #CMAP in rtf
                        continue

                curtag = linee[0].upper()
                continue
            if curtag == "ATOMS":
                if len(linee) >= 4 and linee[0].upper() == "MASS" and isfloat(linee[1]) and isfloat(linee[3]):
                    self.add_atom(linee)
            elif curtag == "BONDS":
                if len(linee) >= 4 and isfloat(linee[2]) and isfloat(linee[3]):
                    self.add_bond(linee)
            elif curtag == "ANGLES":
                if len(linee) >= 5 and isfloat(linee[3]) and isfloat(linee[4]):
                    self.add_angle(linee)
            elif curtag == "DIHEDRALS":
                if len(linee) >= 7 and isfloat(linee[4]) and isfloat(linee[5]) and isfloat(linee[6]):
                    #print linee
                    self.add_dihe(linee)
            elif curtag[0:8] == "IMPROPER":
                if len(linee) >= 7 and isfloat(linee[4]) and isfloat(linee[5]) and isfloat(linee[6]):
                    self.add_impr(linee)
            elif curtag == "NONBONDED":
                if len(linee) >= 4 and isfloat(linee[1]) and isfloat(linee[2]) and isfloat(linee[3]):
                    self.add_nonbond(linee)
            elif curtag == "CMAP":
                self.add_cmap(linee)
                continue
            elif curtag == "NBFIX":
                if len(linee) >= 4 and isfloat(linee[2]) and isfloat(linee[3]):
                    self.add_nbfix(linee)
                continue
            else:
                pass
                #print "skip line", prmfile, line

    def add_atom(self,linee):
        #MASS   192 OC3C61  15.99940 ! ether in six membered ring
        name = linee[2]
        mass = float(linee[3])
        try:
            i = linee.index("!")
            memo = " ".join(linee[i+1:])
        except:
            memo = ""
        atm = prm_atom(name, mass, memo)
        self.atom.append(atm)
        self.atom_dict[name.upper()] = len(self.atom) - 1

    def add_bond(self,linee):
        #                  Kb       k0
        #OC2D3   CC2O3    700.00  1.215   ! ketone MP2/6-31g*, CSD geometry
        atoms = linee[0:2]
        kb = float(linee[2])
        k0 = float(linee[3])
        try:
            i = linee.index("!")
            memo = " ".join(linee[i+1:])
        except:
            memo = ""
        bd = prm_bond(atoms, kb, k0, memo)
        self.bond.append(bd)
        self.bond_dict[" ".join(atoms).upper()] = len(self.bond) - 1
        self.bond_dict[" ".join(atoms[::-1]).upper()] = len(self.bond) - 1

    def add_angle(self,linee):
        #CC312   CC312   CC312    45.00   111.00  ! adm 11/08, glycerol
        #CC3162  CC3161  CC3161   53.35   111.00   8.00  2.561  ! par22 CT1 CT1 CT1
        atoms  = linee[0:3]
        Ktheta = float(linee[3])
        Theta0 = float(linee[4])
        memo = ""
        try:
            Kub = float(linee[5])
            S0  = float(linee[6])
        except:
            Kub = -1
            S0  = -1
        try:
            i = linee.index("!")
            memo = " ".join(linee[i+1:])
        except:
            memo = ""
        ang = prm_angle(atoms, Ktheta, Theta0, Kub, S0, memo)
        self.angle.append(ang)
        name = " ".join(atoms).upper()
        name1 = " ".join(atoms[::-1]).upper()
        if name in self.angle_dict or name1 in self.angle_dict:
            print "WARNING duplicate angle parameter try to add %s %f %f %s"%(name, Ktheta, Theta0, memo)
            if name in self.angle_dict:
                id = self.angle_dict[name]
            if name1 in self.angle_dict:
                id = self.angle_dict[name1]
            print "will replace current value %s %f %f %s"%(" ".join(self.angle[id].atom), \
                    self.angle[id].Ktheta, self.angle[id].Theta0, self.angle[id].memo)
        self.angle_dict[name] = len(self.angle) - 1
        self.angle_dict[name1] = len(self.angle) - 1

    def add_dihe(self,linee):
        #!atom types             Kchi    n   delta
        #CAI  CA   CA   CAI      3.1000  2   180.00 ! from CA CA CA CA
        names = [ i.upper() for i in linee[0:4] ]

        kchi = float(linee[4])
        n = int(linee[5])
        delta = float(linee[6])

        #atom type "X" is replaced with "*"
        #for i in range(0, 4):
        #    if names[i] == "X":
        #        names[i] = "*"

        if delta == 180:
            delta_str = str(0)
        else:
            delta_str = str(int(delta))
        namestr = "_".join(names + [delta_str])
        tag = -1
        for i in range(0, len(self.dihe)):
            if self.dihe[i].namestr == namestr:
                tag = i
                break
        try:
            id = linee.index("!")
            memo = " ".join(linee[id+1:])
        except:
            memo = ""
        if tag == -1:
            dh = prm_dihe(names, kchi, n, delta)
            self.dihe.append(dh)
            dihe_index = len(self.dihe) - 1
            names1 = " ".join(names)
            names1_reverse = " ".join(names[::-1])
            if "X" in names or "*" in names:
                if not names1 in self.dihe_dict["wildcard"]:
                    self.dihe_dict["wildcard"][names1] = []
                if not names1_reverse in self.dihe_dict["wildcard"]:
                    if names1_reverse != names1:
                        self.dihe_dict["wildcard"][names1_reverse] = []
                self.dihe_dict["wildcard"][names1].append(dihe_index)
                if names1_reverse != names1:
                    self.dihe_dict["wildcard"][names1_reverse].append(dihe_index)
            else:
                if not names1 in self.dihe_dict:
                    self.dihe_dict[names1] = []
                if not names1_reverse in self.dihe_dict:
                    if names1_reverse != names1:
                        self.dihe_dict[names1_reverse] = []
                self.dihe_dict[names1].append(dihe_index)
                if names1_reverse != names1:
                    self.dihe_dict[names1_reverse].append(dihe_index)

        else:
            self.dihe[i].add_mult(names, kchi, n, delta)


    def add_impr(self,linee):
        #!atom types           Kpsi                   psi0
        #HR1  NR1  NR2  CPH2    0.5000         0      0.0000 ! ALLOW ARO
        atoms = [n.upper() for n in linee[0:4]]
        #for i in range(0, 4):
        #    if atoms[i].upper() == "X":
        #        atoms[i] = "*"
                
        kpsi = float(linee[4])
        psi0 = float(linee[6])
        try:
            id = linee.index("!")
            memo = " ".join(linee[id+1:])
        except:
            memo = ""

        imp = prm_impr(atoms, kpsi, psi0, memo)
        self.impr.append(imp)
        impindex = len(self.impr) - 1
        names1 = " ".join(atoms)
        names1_reverse = " ".join(atoms[::-1])
        if "X" in atoms or "*" in atoms:
            self.impr_dict["wildcard"][names1] = impindex
        else:
            self.impr_dict[names1] = impindex


    def add_cmap(self,linee):
        try:
            i = linee.index("!")
            memo = " ".join(linee[i+1:])
        except:
            i = 0
            memo = ""

        if len(linee) >= 9 and linee[8].isdigit() == True:
            #C    NH1  CTD1  C    NH1  CTD1  C    NH1   24 
            ntorsion = int(linee[8])
            id = len(self.cmap) + 1
            self.cmap.append(prm_cmap(id, linee[0:8], ntorsion))
            self.cmap_dict[" ".join(linee[0:8]).upper()] = id - 1
            return

        value = []
        if i == 0:
            final = len(linee)
        else:
            final = i
        for j in range(0, final):
            value.append(float(linee[j]))
        self.cmap[-1].add_torsion(value)

    def add_nonbond(self,linee):
        #print linee
        #sigma = exp(log(Rmin) - log(2) / 6)
        #epsilon = - epsilon
        #                  epsilon      Rmin/2                epsilon_14      Rmin/2_14
        #CP1    0.000000  -0.020000     2.275000   0.000000  -0.010000     1.900000 ! ALLOW   ALI
        #CE2    0.000000  -0.064000     2.080000 ! 
        name    = linee[0]
        epsilon = - float(linee[2])
        rmin    = float(linee[3]) * 2
        #['DUM', '0.000000', '-0.000000', '0.000000']
        if rmin > 0:
            sigma   = math.exp(math.log(rmin) - math.log(2)/6)
        else:
            sigma = 0.0
        try:
            epsilon_14 = - float(linee[5])
            rmin_14    = float(linee[6]) * 2
            sigma_14   = math.exp(math.log(rmin_14) - math.log(2)/6)
        except:
            epsilon_14 = -1
            sigma_14   = -1

        try:
            id = linee.index("!")
            memo = " ".join(linee[id+1:])
        except:
            memo = ""

        nb = prm_nonbond(name, epsilon, sigma, epsilon_14, sigma_14, memo)
        self.nonbond.append(nb)
        self.nonbond_dict[name] = len(self.nonbond) - 1
        if name.upper() in ["OT", "HT"]:
            print name, epsilon, rmin, sigma

    def add_nbfix(self, linee):
        #NBFIX atom_i* atom_j*  emin rmin [ emin14 [ rmin14 ]]
        #atom1  atom2   emin rmin
        #SOD    CLA      -0.083875   3.731 !  From osmotic pressure calibration, J. Phys.Chem.Lett. 1:183-189
        #POT    CLA      -0.114236   4.081 !  From osmotic pressure calibration, J. Phys.Chem.Lett. 1:183-189
        name1 = linee[0]
        name2 = linee[1]
        emin = linee[2]
        rmin = linee[3]
        try:
            emin14 = float(linee[4])
        except:
            emin14 = 0
        try:
            rmin14 = float(linee[5])
        except:
            rmin14 = -1
        try:
            id = linee.index("!")
            memo = " ".join(linee[id+1:])
        except:
            memo = ""
        self.nbfix.append(prm_nbfix(name1,name2, emin, rmin, emin14, rmin14, memo))
        name = name1 + " " +name2 
        self.nbfix_dict[name.upper()] = len(self.nbfix) - 1
        name = name2 + " " +name1 
        self.nbfix_dict[name.upper()] = len(self.nbfix) - 1


class dms:
    def __init__(self, dmsfile):
        #self.rules = {}
        #self.rules["info"] = ["prmfile"]
        #self.rules["vdw_func"] = "LJ12_6_sig_epsilon"
        #self.rules["vdw_comb_rule"] = "ARITHMETIC/GEOMETRIC"
        #self.rules["exclusions"] = 4
        #self.rules["es_scale"] = [0.0, 0.0, 1.0]
        #self.rules["lj_scale"] = [0.0, 0.0, 1.0]
        #self.rules["plugins"] = ["bonds", "angles", "ureybradley", "propers", "impropers", "vdw1", "exclusions", "mass"]
        #self.rules["nbfix_identifier"] = "charmm"
        self.dmsfile = dmsfile
        self.conn = sqlite3.connect(self.dmsfile)
        self.c = self.conn.cursor()
        self.c.execute("CREATE TABLE msys_ct ('msys_name' text, id integer primary key)")
        self.c.execute(""" INSERT INTO "msys_ct" VALUES('dms file converted from charmm-gui system. Script written by Yifei Qi qiyifei@gmail.com',0)  """ )

        self.c.execute("CREATE TABLE dms_version (major integer not null, minor integer not null)")
        self.c.execute(""" INSERT INTO "dms_version" VALUES(1 ,7)  """ )
        
        self.c.execute("CREATE TABLE virtual_term (name text);")
        self.c.execute("CREATE TABLE polar_term (name text);")
        self.c.execute("CREATE TABLE provenance(id integer primary key, version text, timestamp text, user text, workdir text, cmdline text, executable text)")
        dttime=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.c.execute(""" INSERT INTO "provenance" VALUES(1, '1.0', '%s', '%s', '%s', '%s', '%s') """%(dttime,getpass.getuser(), os.getcwd(), sys.argv[0], 'charmmgui2dms.py') )

        self.c.execute(""" CREATE TABLE msys_hash('system' text) """)
        hashtag = hash(dttime+sys.argv[0])
        self.c.execute(""" INSERT INTO msys_hash VALUES(%s) """%hashtag)

        self.conn.commit()

    def replace_quota(self, s):
        return s.replace("'", "''")

    def dihe_compare(self, name1, name2):
        match = 0
        for k in range(0, 4):
            if name1[k].upper() == "X" or name2[k].upper() == "X":
                match += 1
                continue
            elif name1[k] == name2[k]:
                match += 1
        if match == 4:
            return True
        #reverse order
        name1.reverse()
        match = 0
        for k in range(0, 4):
            if name1[k].upper() == "X" or name2[k].upper() == "X":
                match += 1
                continue
            elif name1[k] == name2[k]:
                match += 1
        if match == 4:
            return True
        return False

    def generate_dms(self, psf, crd, prm, temperature, boxsize, posref, atom_velocity):
        print "generating",self.dmsfile
        #self.c.execute("""CREATE TABLE forcefield ('path' text, 'info' text)""")
        self.c.execute("""CREATE TABLE forcefield ('path' text)""")
        for ff in prm.prmfiles:
            self.c.execute(""" INSERT INTO "forcefield" VALUES('%s')""" %ff)

        self.c.execute("CREATE TABLE global_cell (id integer primary key,x float, y float, z float)")
        assert len(boxsize) == 3
        count = 0
        for size in boxsize:
            assert len(size) == 3
            count += 1
            self.c.execute("""INSERT INTO "global_cell" VALUES(%d,%f,%f,%f)"""%(count, size[0], size[1], size[2]) )

        self.c.execute("""CREATE TABLE nonbonded_info (vdw_funct text, vdw_rule text,es_funct text)""")
        self.c.execute("""INSERT INTO "nonbonded_info" VALUES('vdw_12_6','arithmetic/geometric','')""")
        
        #nonbond par list
        self.c.execute("CREATE TABLE nonbonded_param ('type' text,\
                'sigma' float, \
                'epsilon' float, \
                'nbfix_identifier' text,\
                'memo' text,\
                id integer primary key) " )
        psfatomtypes = set([i.atomtype for i in psf.atom])
        count = 0
        psfatomtypes_sort = []
        for i in prm.nonbond:
            if i.name in psfatomtypes:
                self.c.execute(""" INSERT INTO "nonbonded_param" VALUES("%s", %.14f, %f, "charmm", "%s", %d) """%\
                        (self.replace_quota(i.name), i.sigma, i.epsilon, self.replace_quota(i.memo), count))
                psfatomtypes_sort.append(i.name)
                count += 1

        #nbfix
        self.c.execute(" CREATE TABLE nonbonded_combined_param ('param1' integer, \
                'param2' integer, \
                'type' text, \
                'sigma' float, \
                'epsilon' float, \
                'nbfix_identifier' text, \
                'memo' text); " )

        for nbfix in prm.nbfix:
            if nbfix.type[0] in psfatomtypes_sort and nbfix.type[1] in psfatomtypes_sort:
                id1 = psfatomtypes_sort.index(nbfix.type[0])
                id2 = psfatomtypes_sort.index(nbfix.type[1])
                names = nbfix.type[0] + " " + nbfix.type[1]
                self.c.execute(""" INSERT INTO "nonbonded_combined_param" VALUES(%d,%d,'%s',%.14f,%f,'charmm','%s') """%\
                        (id1, id2, self.replace_quota(names), nbfix.sigma, nbfix.epsilon_desmond, self.replace_quota(nbfix.memo) ))

        #particle list
        self.c.execute("CREATE TABLE particle (\
            id integer primary key,\
            anum integer,\
            name text not null,\
            x float,\
            y float,\
            z float,\
            vx float,\
            vy float,\
            vz float,\
            resname text not null,\
            resid integer,\
            chain text not null,\
            segid text not null,\
            mass float,\
            charge float,\
            formal_charge integer,\
            resonant_charge float,\
            insertion text not null,\
            msys_ct integer not null,\
            nbtype integer not null\
        )" )
        
        if len(atom_velocity) == 0:
            #velocity generation; standard normal distribution
            atom_velocity = np.random.normal(loc=0.0, scale=1, size = (psf.atomnum, 3)) 
            atommass = np.array([i.mass for i in psf.atom])
            vstd = 151.3343 * np.sqrt(temperature/atommass) * 0.01
            #change std
            atom_velocity = np.multiply(atom_velocity.T, vstd).T

            #remove COM move
            momento = np.multiply(atom_velocity.T, atommass).T
            mean_velo = np.sum(momento, axis = 0)/ np.sum(atommass)
            atom_velocity = atom_velocity - mean_velo

        for i in range(0, psf.atomnum):
            primary_key = i
            anum = psf.atom[i].anum
            name = psf.atom[i].atomname
            #to escape "'" in atom name
            #name = name.replace("'", "''")
            mass = psf.atom[i].mass
            x, y, z = crd.get_coor(i)
            vx, vy, vz = atom_velocity[i]

            resname = psf.atom[i].resname
            resid = psf.atom[i].resid
            chainid = psf.atom[i].segid
            segid = psf.atom[i].segid
            charge = psf.atom[i].charge
            formal_charge = 0
            resonant_charge = 0 
            if resid.isdigit() == False:
                #resid has insertion, e.g. 112A
                insertion = resid[-1]
                resid = int(resid[0:-1])
            else:
                resid = int(resid)
                insertion = ""

            msys_ct = 0
            nbtype = psfatomtypes_sort.index(psf.atom[i].atomtype)
            self.c.execute(""" INSERT INTO "particle" VALUES(%d,%d,'%s',%f,%f,%f,%.14f,%.14f,%.14f,'%s',%d,'%s','%s',%f,%f,%f,%f,'%s',%d,%d)"""%(i, anum, self.replace_quota(name), x, y, z, vx, vy, vz, self.replace_quota(resname), resid, self.replace_quota(chainid), self.replace_quota(segid), mass, charge, formal_charge, resonant_charge, insertion, msys_ct, nbtype) )
        print "particle table done"

        #bond
        self.c.execute("CREATE TABLE bond_term (name text)")
        if len(psf.bond) > 0:
            self.c.execute(""" INSERT INTO "bond_term" VALUES('stretch_harm') """ )

        #TIP3-H -- TIP3-H bond is listed in psf file but is not included in psf.bond
        self.c.execute("CREATE TABLE bond ( p0 integer, p1 integer, 'order' interger, resonant_order float)")
        for bond in psf.bond:
            id1 = bond[0] - 1
            id2 = bond[1] - 1
            if id1 < id2:
                self.c.execute(""" INSERT INTO "bond" VALUES(%d, %d, %d, %f) """ %(id1, id2, 1.0, 1.0))
            else:
                self.c.execute(""" INSERT INTO "bond" VALUES(%d, %d, %d, %f) """ %(id2, id1, 1.0, 1.0))

        self.c.execute("CREATE TABLE stretch_harm_param ('type' text, 'r0' float, 'fc' float, 'memo' text, id integer primary key)")
        bond_par_id = {}
        for i in range(0, len(prm.bond)):
            type = prm.bond[i].atom[0] + " " + prm.bond[i].atom[1]
            type1 = prm.bond[i].atom[1] + " " + prm.bond[i].atom[0]
            bond_par_id[type] = i
            bond_par_id[type1] = i
            memo = prm.bond[i].memo.replace("'","")
            self.c.execute(""" INSERT INTO "stretch_harm_param" VALUES('%s',%f,%f,'%s',%d) """%(self.replace_quota(type), prm.bond[i].b0, prm.bond[i].kb, self.replace_quota(memo), i) )
        bond_prm_number = i 
        
        self.c.execute("CREATE TABLE stretch_harm_term (p0 integer, p1 integer, 'constrained' integer, param integer not null)")
        for bond in psf.bond:
            id1 = bond[0] - 1
            id2 = bond[1] - 1
            type1 = psf.atom[id1].atomtype
            type2 = psf.atom[id2].atomtype
            type = type1 + " " + type2
            parid = bond_par_id[type]
            if psf.atom[id1].anum == 1 or psf.atom[id2].anum == 1: #bonds having hydrogen are constrained
                constraint = 1
            else:
                constraint = 0
            if id1 < id2:
                self.c.execute(""" INSERT INTO "stretch_harm_term" VALUES(%d,%d,%d,%d) """%(id1, id2, constraint, parid) )
            else:
                self.c.execute(""" INSERT INTO "stretch_harm_term" VALUES(%d,%d,%d,%d) """%(id2, id1, constraint, parid) )
            #UB terms are added in angle section

        print "bond table done"
        #angles
        if len(psf.angle) > 0:
            self.c.execute(""" INSERT INTO "bond_term" VALUES('angle_harm')""" )
        self.c.execute("CREATE TABLE angle_harm_term (p0 integer, p1 integer, p2 integer, 'constrained' integer, param integer not null)" )
        self.c.execute("CREATE TABLE angle_harm_param ('type' text, 'theta0' float, 'fc' float, 'memo' text, id integer primary key)" )

        psf_angles = set()
        for angle in psf.angle:
            id1 = angle[0]
            id2 = angle[1]
            id3 = angle[2]
            type1 = psf.atom[id1-1].atomtype
            type2 = psf.atom[id2-1].atomtype
            type3 = psf.atom[id3-1].atomtype
            pair1 = " ".join([type1, type2, type3])
            pair2 = " ".join([type3, type2, type1])
            psf_angles.add(pair1)
            psf_angles.add(pair2)
        angle_prm_id = 0
        angle_pair2id = {}
        ub_pair2id = {}
        for angle in prm.angle:
            pair = " ".join(angle.atom).upper()
            pair1 = " ".join(angle.atom[::-1]).upper()
            if pair in psf_angles or pair1 in psf_angles:
                angle_pair2id[pair] = angle_prm_id
                angle_pair2id[pair1] = angle_prm_id
                self.c.execute(""" INSERT INTO "angle_harm_param" VALUES('%s',%f,%f,'%s',%d) """ %(self.replace_quota(pair), angle.Theta0, angle.Ktheta, self.replace_quota(angle.memo), angle_prm_id))
                angle_prm_id += 1
                #UB term
                if angle.Kub != -1 and angle.S0 != -1:
                    bond_prm_number += 1
                    type = "ureybradley " + pair
                    type1 = "ureybradley " + pair1
                    self.c.execute(""" INSERT INTO "stretch_harm_param" VALUES('%s',%f,%f,'%s',%d) """%(self.replace_quota(type), angle.S0, angle.Kub,"ureybradley "+pair+" "+ self.replace_quota(angle.memo), bond_prm_number) )
                    ub_pair2id[type] = bond_prm_number
                    ub_pair2id[type1] = bond_prm_number

        for angle in psf.angle:
            id1 = angle[0]
            id2 = angle[1]
            id3 = angle[2]
            type1 = psf.atom[id1-1].atomtype
            type2 = psf.atom[id2-1].atomtype
            type3 = psf.atom[id3-1].atomtype
            pair = " ".join([type1, type2, type3]).upper()
            prmid = angle_pair2id[pair]
            #TIP3 angle is constrained
            if psf.atom[id1-1].resname == "TIP3" and psf.atom[id2-1].resname == "TIP3" and psf.atom[id2-1].resname == "TIP3":
                cons = 1
            else:
                cons = 0

            self.c.execute(""" INSERT INTO "angle_harm_term" VALUES(%d,%d,%d,%d,%d) """%(id1-1, id2-1, id3-1, cons, prmid))
            ubtype = "ureybradley "+pair
            #this angle has UB term
            if ubtype in ub_pair2id:
                constraint = 0 #is constrained?? TBD
                self.c.execute(""" INSERT INTO "stretch_harm_term" VALUES(%d,%d,%d,%d) """%(id1-1, id3-1, constraint, ub_pair2id[ubtype]) )

        print "angle table done"
        #dihedrals   
        if len(psf.dihe) > 0:
            self.c.execute(""" INSERT INTO "bond_term" VALUES('dihedral_trig') """ )
        dihpair_psf = psf.get_uniq_dihedral_pairs()
        dihpair_psf1 = [i.split() for i in dihpair_psf]

        self.c.execute(""" CREATE TABLE dihedral_trig_param ('type' text, 'phi0' float, 'fc0' float, 'fc1' float, 'fc2' float, 'fc3' float, 'fc4' float, 'fc5' float, 'fc6' float, 'memo' text, id integer primary key)""" )

        prm_id = -1
        dihe_prm_uniq = []
        dihe2id = {}
        for i in prm.dihe:
            for dihe in dihpair_psf1:
                match = self.dihe_compare(i.names, dihe)
                if match:
                    #print i.names, i.delta, i.fc, "matched to",dihe
                    dihe_prm_type = " ".join(i.names).upper()
                    if not (dihe_prm_type, i.delta) in dihe_prm_uniq:
                        dihe_prm_uniq.append((dihe_prm_type, i.delta))
                        prm_id += 1
                        self.c.execute(""" INSERT INTO "dihedral_trig_param" VALUES('%s',%f,%f,%f,%f,%f,%f,%f,%f,'%s',%d) """%(self.replace_quota(dihe_prm_type), i.delta,\
                                i.fc[0], i.fc[1], i.fc[2], i.fc[3], i.fc[4], i.fc[5], i.fc[6], self.replace_quota(i.memo), prm_id))
                    if not " ".join(dihe).upper() in dihe2id:
                        dihe2id[" ".join(dihe).upper()] = [prm_id]
                    else:
                        dihe2id[" ".join(dihe).upper()].append(prm_id)

        self.c.execute(""" CREATE TABLE dihedral_trig_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null) """ )
        for dihe in psf.dihe:
            types = [ psf.atom[j-1].atomtype.upper() for j in dihe ]
            type1 = " ".join(types)
            if not type1 in dihe2id:
                print "can not find parameter for this dihedral"
                for j in dihe:
                    print j, psf.atom[j-1].atomstr

            prmids = dihe2id[type1]
            prmids1 = []
            if len(prmids) == 1:
                prmids1 = prmids
            else:
                #if more than 1 match, remove wildcard match
                for j in prmids:
                    if not ("X" in dihe_prm_uniq[j][0].split() or "*" in dihe_prm_uniq[j][0].split()):
                        prmids1.append(j)
            assert len(prmids1) > 0

            for prmid in prmids1:
                self.c.execute(""" INSERT INTO "dihedral_trig_term" VALUES(%d,%d,%d,%d,%d) """%(dihe[0]-1, dihe[1]-1, dihe[2]-1, dihe[3]-1, prmid ))

        print "dihedral table done"
        #improper
        if len(psf.impr) > 0:
            self.c.execute(""" CREATE TABLE improper_harm_param ('type' text, 'phi0' float, 'fc' float, 'memo' text, id integer primary key) """ )
            self.c.execute(""" INSERT INTO "bond_term" VALUES('improper_harm') """ )
            self.c.execute(""" CREATE TABLE improper_harm_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null)""" )

            psf_imprpair = psf.get_uniq_impr_pairs()
            prm_psf_impr = []
            impr2id = {}
            for prm_impr in prm.impr:
                type = " ".join(prm_impr.atoms)
                for psf_impr in psf_imprpair:
                    match = self.dihe_compare(prm_impr.atoms, psf_impr.split())
                    if match:
                        prm_psf_impr.append(type)
                        id = len(prm_psf_impr) - 1
                        if not psf_impr in impr2id:
                            impr2id[psf_impr] = [id]
                        else:
                            impr2id[psf_impr].append(id)
                        self.c.execute(""" INSERT INTO "improper_harm_param" VALUES('%s',%f,%f,'%s',%d) """ %(self.replace_quota(type), prm_impr.psi0\
                            , prm_impr.kpsi, prm_impr.memo, id))

            for psf_impr in psf.impr:
                type = " ".join([psf.atom[i-1].atomtype for i in psf_impr])
                if not type in impr2id:
                    print psf_impr, type
                ids = impr2id[type]
                ids1 = []
                if len(ids) == 1:
                    ids1 = ids
                else:
                    #if more than 1 match, remove wildcard match
                    for j in ids:
                        if prm_psf_impr[j].count(" X ") == 0:
                            ids1.append(j)
                assert len(ids1) > 0
                if len(ids1) > 1:
                    print type, "has more than 1 impr"
                for id in ids1:
                    self.c.execute(""" INSERT INTO "improper_harm_term" VALUES(%d,%d,%d,%d,%d)""" %(psf_impr[0]-1, \
                            psf_impr[1]-1, psf_impr[2]-1, psf_impr[3]-1, id) )

        print "improper table done"
        #CMAP
        if len(psf.cmap) > 0:
            self.c.execute(""" INSERT INTO "bond_term" VALUES('torsiontorsion_cmap') """ )
            self.c.execute(""" CREATE TABLE torsiontorsion_cmap_term (p0 integer, p1 integer, p2 integer, p3 integer, p4 integer, p5 integer, p6 integer, p7 integer, param integer not null)""" )
            self.c.execute(""" CREATE TABLE torsiontorsion_cmap_param ('type' text, 'cmapid' text, 'memo' text, id integer primary key)""" )

            id = 0
            cmap2id = {}
            for cmap in prm.cmap:
                self.c.execute(""" INSERT INTO "torsiontorsion_cmap_param" VALUES('%s','cmap%d','%s',%d)"""%(self.replace_quota(" ".join(cmap.atom)), cmap.id, "", id) )
                cmap2id[" ".join(cmap.atom).upper()] = id
                self.c.execute(""" CREATE TABLE cmap%d ('phi' float, 'psi' float, 'energy' float)""" %cmap.id)
                energyid = 0
                for phi in cmap.phi:
                    for psi in cmap.psi:
                        self.c.execute(""" INSERT INTO "cmap%d" VALUES(%f,%f,%f)"""%(cmap.id, phi, psi, cmap.cmap[energyid]) )
                        energyid += 1
                id += 1

            for cmap in psf.cmap:
                atomnames = [psf.atom[i-1].atomtype.upper() for i in cmap]
                cmap_type = " ".join(atomnames) 
                cmapid = cmap2id[cmap_type]
                self.c.execute(""" INSERT INTO "torsiontorsion_cmap_term" VALUES(%d,%d,%d,%d,%d,%d,%d,%d,%d)""" %(cmap[0]-1, \
                        cmap[1]-1, cmap[2]-1, cmap[3]-1, cmap[4]-1, cmap[5]-1, cmap[6]-1, cmap[7]-1, cmapid ))

        print "CMAP table done"
        #1-4 pairs
        n14pair = 0
        if len(psf.dihe) > 0:
            #get 1-4 pairs from dihedral list
            #do not include if the pair is in a ring and also forms a bond or angle, e.g., CD2 and NE1 in TRP
            self.c.execute(""" INSERT INTO "bond_term" VALUES('pair_12_6_es')""" )
            self.c.execute(""" CREATE TABLE pair_12_6_es_term (p0 integer, p1 integer, param integer not null)""" )

            self.c.execute("SELECt p0, p1 from stretch_harm_term")
            bondpair = self.c.fetchall()
            self.c.execute("SELECt p0, p2 from angle_harm_term")
            anglepair = self.c.fetchall()
            no14pair = set(bondpair + anglepair)
                    
            self.c.execute("SELECt p0, p3 from dihedral_trig_term")
            dihepair = set(self.c.fetchall())
            uniq_14pair = []
            type2prmid = {}
            for dihe in dihepair:
                id1 = dihe[0]
                id2 = dihe[1]
                type1 = psf.atom[id1].atomtype
                type2 = psf.atom[id2].atomtype
                charge1 = psf.atom[id1].charge
                charge2 = psf.atom[id2].charge
                if (id1, id2) in no14pair or (id2, id1) in no14pair:
                    continue
                else:
                    if not ( (type1, charge1, type2, charge2) in uniq_14pair or (type2, charge2, type1, charge1) in uniq_14pair):
                        uniq_14pair.append((type1, charge1, type2, charge2))
                        type2prmid[(type1, charge1, type2, charge2)] = len(uniq_14pair) - 1
                        type2prmid[(type2, charge2, type1, charge1)] = len(uniq_14pair) - 1

                    parid=type2prmid[(type1, charge1, type2, charge2)]
                    self.c.execute(""" INSERT INTO "pair_12_6_es_term" VALUES(%d,%d,%d)"""%(id1, id2, parid) )
                    n14pair += 1
            print "found %d 1-4 pairs"%n14pair
            
            self.c.execute(""" CREATE TABLE pair_12_6_es_param ('aij' float, 'bij' float, 'qij' float, 'type' text, 'memo' text, id integer primary key)""")

            for i in range(0, len(uniq_14pair)):
                type1 = uniq_14pair[i][0]
                charge1 = uniq_14pair[i][1]
                type2 = uniq_14pair[i][2]
                charge2 = uniq_14pair[i][3]
                e1 = prm.nonbond[prm.nonbond_dict[type1]].epsilon_14
                s1 = prm.nonbond[prm.nonbond_dict[type1]].sigma_14
                if s1 < 0:
                    #print type1, "does not have 14 interaction para, copy from nonbond para"
                    e1 =  prm.nonbond[prm.nonbond_dict[type1]].epsilon
                    s1 = prm.nonbond[prm.nonbond_dict[type1]].sigma
                e2 = prm.nonbond[prm.nonbond_dict[type2]].epsilon_14
                s2 = prm.nonbond[prm.nonbond_dict[type2]].sigma_14
                if s2 < 0:
                    #print type2, "does not have 14 interaction para, copy from nonbond para"
                    e2 =  prm.nonbond[prm.nonbond_dict[type2]].epsilon
                    s2 = prm.nonbond[prm.nonbond_dict[type2]].sigma

                #arithmetic average on sigma , and geometric average on epsilon
                sigma = (s1 + s2) / 2
                epsilon = math.sqrt(e1 * e2)
                #desmond manual 7.3.1
                aij = 4*epsilon*sigma**12
                #bij is positive, verified
                bij = 4*epsilon*sigma**6
                assert aij > 0
                qij = charge1 * charge2
                self.c.execute(""" INSERT INTO "pair_12_6_es_param" VALUES(%.14f,%.14f,%f,'%s','%s',%d)"""%(aij,  bij, qij, self.replace_quota(type1+" "+type2+" 14 inter"), "", i) )
        print "1-4 pair table done"

        #exclusion list
        self.c.execute("select p0, p1 from bond")
        bondpair = self.c.fetchall()
        self.c.execute("select p0, p2 from angle_harm_term")
        anglepair = self.c.fetchall()
        self.c.execute("select p0, p3 from dihedral_trig_term")
        dihepair = self.c.fetchall()
        exclusion_pair = set()
        for i in bondpair + anglepair + dihepair:
            if i[0] > i[1]:
                exclusion_pair.add((i[1], i[0]))
            else:
                exclusion_pair.add((i[0], i[1]))
        print "found",len(exclusion_pair), "exclusion pairs"
        if len(exclusion_pair) > 0:
            self.c.execute("CREATE TABLE exclusion (p0 integer, p1 integer)")
            for i in exclusion_pair:
                self.c.execute("""INSERT INTO "exclusion" VALUES(%d,%d)"""%(i[0], i[1]))

        print "exclusion table done"

        #constraints
        self.c.execute("CREATE TABLE constraint_term (name text)")
        #one hydrogen atom
        self.c.execute("""INSERT INTO "constraint_term" VALUES('constraint_ah1')""")
        self.c.execute("""CREATE TABLE constraint_ah1_term (p0 integer, p1 integer, param integer not null)""")
        self.c.execute("""CREATE TABLE constraint_ah1_param ('r1' float, id integer primary key)""")

        #two hydrogen atom
        self.c.execute("""INSERT INTO "constraint_term" VALUES('constraint_ah2')""")
        self.c.execute("""CREATE TABLE constraint_ah2_term (p0 integer, p1 integer, p2 integer, param integer not null)""")
        self.c.execute("""CREATE TABLE constraint_ah2_param ('r1' float, 'r2' float, id integer primary key)""")

        #three hydrogen atom
        self.c.execute("""INSERT INTO "constraint_term" VALUES('constraint_ah3')""")
        self.c.execute("""CREATE TABLE constraint_ah3_term (p0 integer, p1 integer, p2 integer, p3 integer, param integer not null)""")
        self.c.execute("""CREATE TABLE constraint_ah3_param ('r1' float, 'r2' float, 'r3' float, id integer primary key)""")

        #TIP3
        self.c.execute("""INSERT INTO "constraint_term" VALUES('constraint_hoh')""")
        self.c.execute("""CREATE TABLE constraint_hoh_term (p0 integer, p1 integer, p2 integer, param integer not null)""")
        self.c.execute("""CREATE TABLE constraint_hoh_param ('theta' float, 'r1' float, 'r2' float, id integer primary key)""")

        #count how many hydrogen atoms each heavy atom is bonded
        nhydrogen_dict = {}
        self.c.execute("select id from particle where anum > 1")
        for i in self.c.fetchall():
            #this dict containts the id of hydrogens that are bonded 
            nhydrogen_dict[i[0]] = []

        for bondpair in psf.bond:
            anum1 = psf.atom[bondpair[0]-1].anum
            anum2 = psf.atom[bondpair[1]-1].anum
            if anum1 != 1 and anum2 == 1:
                nhydrogen_dict[bondpair[0]-1].append(bondpair[1]-1)
            elif anum1 == 1 and anum2 != 1:
                nhydrogen_dict[bondpair[1]-1].append(bondpair[0]-1)

        constraint_para_dict = {}
        for i in range(1, 5):
            #store uniq para set for each Ahn constarint type
            constraint_para_dict[i] = []
        constraint_para_dict[-1] = [] #for water parameter

        for heavyatom_id in nhydrogen_dict:
            nhydrogen = len(nhydrogen_dict[heavyatom_id]) 
            if nhydrogen == 0:
                continue
            if nhydrogen > 4:
                raise Exception("more than 4 hydrogens bonded to one heavy atom %s is not support"%psf.atom[heavyatom_type].atomstr)
            heavyatom_type = psf.atom[heavyatom_id].atomtype
            hydrogenatom_type = [psf.atom[i].atomtype for i in nhydrogen_dict[heavyatom_id]]
            bondpar_index = [ prm.get_bond_par_index(heavyatom_type,i) for i in hydrogenatom_type]
            bondlength = [prm.bond[i].b0 for i in bondpar_index]

            if psf.atom[heavyatom_id].resname.upper() == "TIP3" and nhydrogen == 2:
                angle_par_index = prm.get_angle_par_index(hydrogenatom_type[0], heavyatom_type, hydrogenatom_type[1])
                parset = [prm.angle[angle_par_index].Theta0, bondlength[0], bondlength[1]]
                if not parset in constraint_para_dict[-1]:
                    constraint_para_dict[-1].append(parset)
                parid = constraint_para_dict[-1].index(parset)
                self.c.execute("""insert into "constraint_hoh_term" values(%d, %d, %d, %d)"""%(heavyatom_id,\
                        nhydrogen_dict[heavyatom_id][0], nhydrogen_dict[heavyatom_id][1], parid))
            else:
                if not bondlength in constraint_para_dict[nhydrogen]:
                    constraint_para_dict[nhydrogen].append(bondlength)
                parid = constraint_para_dict[nhydrogen].index(bondlength)
                
                atomids = [str(heavyatom_id)]
                atomids = atomids + [str(i) for i in nhydrogen_dict[heavyatom_id]]
                atomids = ",".join(atomids)
                self.c.execute("""INSERT into "constraint_ah%d_term" values(%s, %d)"""%(nhydrogen,atomids, parid))

        for n in constraint_para_dict:
            count = 0
            for par in constraint_para_dict[n]:
                if n == -1:
                    self.c.execute("""insert into "constraint_hoh_param" values(%f, %f, %f, %d)"""%(par[0], par[1], par[2], count))
                else:
                    parstr = [str(i) for i in par]
                    parstr = ",".join(parstr)
                    self.c.execute("""insert into "constraint_ah%d_param" values(%s, %d)"""%(n, parstr, count))
                count += 1

        print "constraints done"
        #views
        self.c.execute("""CREATE VIEW angle_harm as 
          select p0, p1, p2, 
        "type", "theta0", "fc", "memo", "constrained"  from angle_harm_param
          join angle_harm_term
          on param=id;""")
        self.c.execute("""CREATE VIEW constraint_ah1 as 
          select p0, p1, 
        "r1"  from constraint_ah1_param
          join constraint_ah1_term
          on param=id;""")
        self.c.execute("""CREATE VIEW constraint_ah2 as 
          select p0, p1, p2, 
        "r1", "r2"  from constraint_ah2_param
          join constraint_ah2_term
          on param=id;""")
        self.c.execute("""CREATE VIEW constraint_ah3 as 
          select p0, p1, p2, p3, 
        "r1", "r2", "r3"  from constraint_ah3_param
          join constraint_ah3_term
          on param=id;""")
        self.c.execute("""CREATE VIEW constraint_hoh as 
          select p0, p1, p2, 
        "theta", "r1", "r2"  from constraint_hoh_param
          join constraint_hoh_term
          on param=id;""")
        self.c.execute("""CREATE VIEW dihedral_trig as 
          select p0, p1, p2, p3, 
        "type", "phi0", "fc0", "fc1", "fc2", "fc3", "fc4", "fc5", "fc6", "memo"  from dihedral_trig_param
          join dihedral_trig_term
          on param=id;""")
        if len(psf.impr) > 0:
            self.c.execute("""CREATE VIEW improper_harm as 
            select p0, p1, p2, p3, 
            "type", "phi0", "fc", "memo"  from improper_harm_param
            join improper_harm_term
            on param=id;""")
        if n14pair > 0:
            self.c.execute("""CREATE VIEW pair_12_6_es as 
            select p0, p1, 
            "aij", "bij", "qij", "type", "memo"  from pair_12_6_es_param
            join pair_12_6_es_term
            on param=id;""")
        self.c.execute("""CREATE VIEW stretch_harm as 
          select p0, p1, 
        "type", "r0", "fc", "memo", "constrained"  from stretch_harm_param
          join stretch_harm_term
          on param=id;""")
        if len(psf.cmap) > 0:
            self.c.execute("""CREATE VIEW torsiontorsion_cmap as 
            select p0, p1, p2, p3, p4, p5, p6, p7, 
            "type", "cmapid", "memo"  from torsiontorsion_cmap_param
            join torsiontorsion_cmap_term
            on param=id;""")


        if len(posref) > 0:
            print "generating position restraints"
            self.c.execute(""" INSERT INTO "bond_term" VALUES('posre_harm') """ )
            self.c.execute("""CREATE TABLE posre_harm_term (p0 integer, 'x0' float, 'y0' float, 'z0' float, param integer not null)""")
            self.c.execute("""CREATE TABLE posre_harm_param ('fcx' float, 'fcy' float, 'fcz' float, id integer primary key)""")
            parid = 0
            for ref in posref:
                force = ref.posk
                self.c.execute(""" insert into "posre_harm_param" values(%f, %f, %f, %d) """%(force[0], force[1], force[2], parid) )
                for atom in ref.atom:
                    id = atom.atomid - 1
                    xyz = atom.coor
                    self.c.execute("""insert into "posre_harm_term" values(%d, %f, %f, %f, %d) """%(id, xyz[0], xyz[1], xyz[2], parid))
                parid += 1

            self.c.execute(""" CREATE VIEW posre_harm as 
              select p0, "fcx", "fcy", "fcz", "x0", "y0", "z0"  from posre_harm_param
              join posre_harm_term on param=id; """)

        self.conn.commit()

class namdvel:
    def __init__(self, vfile):
        print vfile
        f = file(vfile, mode='r')
        self.natom = int(np.fromfile(f, dtype='i4', count=1))
        count = 0
        pdbvelfactor=20.45482706
        self.velocity = np.fromfile(f, dtype='f8', count=3 * self.natom) * pdbvelfactor
        self.velocity = self.velocity.reshape((self.natom, 3))
        print "read velocity for %d atoms from %s"%(self.natom, vfile)

class namdcoor:
    def __init__(self, cfile):
        print cfile
        f = file(cfile, mode='r')
        self.natom = int(np.fromfile(f, dtype='i4', count=1))
        self.coor = np.fromfile(f, dtype='f8', count=3 * self.natom)
        self.coor = self.coor.reshape((self.natom, 3))
        print "read coordinates for %d atoms from %s"%(self.natom, cfile)

    def get_coor(self, i):
        return self.coor[i,:]

class lammps:
    def __init__(self, filename, title="  "):
        self.fp = open(filename, "w")
        self.fp.write("%s\n\n"%title)

    def write(self, line):
        self.fp.write("%s\n"%line)

    def generate_lammps(self, psf, crd, prm, boxsize, lammpscmap, lammpsangle, zcenter):
        #do the dihedrals first because we need an accurate count of dihedral types

        #dihedrals
        # k n d weight
        #determine the weight of each dihe prm
        #if the 1st and 4th atoms form a bond or angle, weight is 0
        #if the 1st and 4th atoms form another dihedral, weight is 0.5
        anglepair = [ (i[0], i[2]) for i in psf.angle]
        bondpair = [ (i[0], i[1]) for i in psf.bond]
        anglebondpair = anglepair + bondpair
        anglebondpair = set(anglebondpair)

        dihe_14atom_dict = {}
        for i in range(0, psf.dihenum):
            id1 = psf.dihe[i][0]
            id2 = psf.dihe[i][3]
            key = (id1, id2)
            key1 = (id2, id1)
            if not key in dihe_14atom_dict:
                dihe_14atom_dict[key] = []
            if not key1 in dihe_14atom_dict:
                dihe_14atom_dict[key1] = []
            dihe_14atom_dict[key].append(i)
            dihe_14atom_dict[key1].append(i)


        if psf.dihenum > 0:
            #"Dihedrals" lines
            lammps_dihe_lines = []
            #"Dihedrals Coeffs" lines
            lammps_dihe_pars = []
            lammps_dihe_id = {}
            for i in range(0, psf.dihenum):
                id1 = psf.dihe[i][0]
                id2 = psf.dihe[i][3]
                atomnames = [psf.atom[j-1].atomname+"-"+psf.atom[j-1].atomtype for j in psf.dihe[i]]
                atomnames = " ".join(atomnames)
                if (id1, id2) in anglebondpair or (id2, id1) in anglebondpair:
                    weight = 0
                elif len(dihe_14atom_dict[(id1, id2)]) == 2:
                    weight = 0.5
                elif len(dihe_14atom_dict[(id1, id2)]) > 2:
                    raise Exception("atom %d and %d apprear in 1-4 atoms in more than three dihedrals"%(id1, id2))
                else:
                    weight = 1

                for prmid in psf.dihe_parid[i]:
                    if psf.dihe_parid[i].index(prmid) > 0:
                        #if matched to mutiple dihe pars, only the first one has weight 1
                        weight = 0
                    if (prmid, weight) in lammps_dihe_id:
                        #This is the ids in the Dihedral Coeffs in the output file
                        ids = lammps_dihe_id[(prmid, weight)]
                    else:
                        prmlines = prm.dihe[prmid].generate_lammps_coef(weight)
                        n1 = len(lammps_dihe_pars)
                        lammps_dihe_pars.extend(prmlines)
                        n2 = len(lammps_dihe_pars)
                        ids = range(n1+1, n2+1)
                        lammps_dihe_id[(prmid, weight)] = ids

                    for id in ids:
                        #this is the dihedral line in the output file
                        lammps_dihe_lines.append("%8d %12d %12d %12d %12d # %s"%(id, psf.dihe[i][0], psf.dihe[i][1], psf.dihe[i][2], psf.dihe[i][3], atomnames))
        
        if psf.atomnum > 0:
            self.write("%18d  atoms"%psf.atomnum)
        if psf.bondnum > 0:
            self.write("%18d  bonds"%psf.bondnum)
        if psf.anglenum > 0:
            self.write("%18d  angles"%psf.anglenum)
        if psf.dihenum > 0:
            self.write("%18d  dihedrals"%len(lammps_dihe_lines))
        if psf.imprnum > 0:
            self.write("%18d  impropers"%psf.imprnum)
        if psf.cmapnum > 0:
            self.write("%18d  crossterms"%psf.cmapnum)
        self.write("\n")

        nonbondtypes = set(psf.nonbond_parid)
        bondtypes = set(psf.bond_parid)
        angletypes = set(psf.angle_parid)
        imprtypes = set(psf.impr_parid)

        self.write("%18d  atom types"%len(nonbondtypes))
        if len(bondtypes) > 0:
            self.write("%18d  bond types"%len(bondtypes))
        if len(angletypes) > 0:
            self.write("%18d  angle types"%len(angletypes))
        if psf.dihenum > 0:
            self.write("%18d  dihedral types"%(len(lammps_dihe_pars)))
        if len(imprtypes) > 0:
            self.write("%18d  improper types\n"%len(imprtypes))

        self.write("%12.4f %12.4f xlo xhi"%(-boxsize[0][0]/2.0, boxsize[0][0]/2.0))
        self.write("%12.4f %12.4f ylo yhi"%(-boxsize[1][1]/2.0, boxsize[1][1]/2.0))
        self.write("%12.4f %12.4f zlo zhi\n"%(-boxsize[2][2]/2.0 + zcenter, boxsize[2][2]/2.0 + zcenter))

        self.write("Masses\n")

        nonbondtypes = list(nonbondtypes)
        nonbondtypes.sort()
        nonbondatomtypes = []
        for j in range(0, len(nonbondtypes)):
            i = nonbondtypes[j]
            atomtype = prm.nonbond[i].name
            nonbondatomtypes.append(atomtype.upper())
            massid = prm.get_mass_par_index(atomtype)
            mass = prm.atom[massid].mass
            self.write("%8d %13.3f  # %s"%(j+1, mass, atomtype))
        self.write("\n")
        
       
        #if the sytem has NBfix, we use PairIJ Coeffs to list all possible pairs
        nbfixtype = []
        for i in range(0, len(prm.nbfix)):
            type1 = prm.nbfix[i].type[0].upper()
            type2 = prm.nbfix[i].type[1].upper()
            if type1 in nonbondatomtypes and type2 in nonbondatomtypes:
                nbfixtype.append([type1, type2])
                nbfixtype.append([type2, type1])

        count = 1
        #have NBfix
        if len(nbfixtype) > 0:
            self.write("\nPairIJ Coeffs\n")
            #nonbondtypes list the uniq list index in prm.nonbond
            for ii in range(0, len(nonbondtypes)):
                for jj in range(ii, len(nonbondtypes)):
                    i = nonbondtypes[ii]
                    j = nonbondtypes[jj]
                    type1 = prm.nonbond[i].name.upper()
                    type2 = prm.nonbond[j].name.upper()
                    type = type1 + " " + type2
                    if [type1, type2] in nbfixtype:
                        type = type + " NBfix"
                        nbfixid = prm.get_nbfix_par_index(type1, type2)
                        if nbfixid < 0:
                            raise Exception("look for nbfix type error %s %s"%(type1, type2))
                        sigma = prm.nbfix[nbfixid].sigma
                        epsilon = math.fabs(prm.nbfix[nbfixid].emin)
                        if prm.nbfix[nbfixid].sigma14 < 0:
                            sigma14 = sigma
                            epsilon14 = epsilon
                        else:
                            sigma14 = prm.nbfix[nbfixid].sigma14
                            epsilon14 = math.fabs(prm.nbfix[nbfixid].emin14)
                    else:
                        sigma = (prm.nonbond[i].sigma + prm.nonbond[j].sigma)*0.5
                        epsilon = math.sqrt(prm.nonbond[i].epsilon*prm.nonbond[j].epsilon)
                        if prm.nonbond[i].sigma_14 < 0:
                            sigma14a = prm.nonbond[i].sigma
                            epsilon14a = prm.nonbond[i].epsilon
                        else:
                            sigma14a = prm.nonbond[i].sigma_14
                            epsilon14a = prm.nonbond[i].epsilon_14
                        if prm.nonbond[j].sigma_14 < 0:
                            sigma14b = prm.nonbond[j].sigma
                            epsilon14b = prm.nonbond[j].epsilon
                        else:
                            sigma14b = prm.nonbond[j].sigma_14
                            epsilon14b = prm.nonbond[j].epsilon_14
                        sigma14 = (sigma14a + sigma14b) * 0.5
                        epsilon14 = math.sqrt(epsilon14a * epsilon14b)

                    self.write("%8d %8d %.14f %.14f %.14f %.14f # %s"%(ii+1, jj+1, epsilon, sigma, epsilon14, sigma14, type))
                    count += 1
        else:
            #if not NBfix, we use Pair Coeffs
            self.write("Pair Coeffs\n")
            for j in range(0, len(nonbondtypes)):
                i = nonbondtypes[j]
                if prm.nonbond[i].sigma_14 < 0:
                    sigma14 = prm.nonbond[i].sigma
                    epsilon14 = prm.nonbond[i].epsilon
                else:
                    sigma14 = prm.nonbond[i].sigma_14
                    epsilon14 = prm.nonbond[i].epsilon_14

                self.write("%8d %.14f %.14f %.14f %.14f # %s"%(j+1, math.fabs(prm.nonbond[i].epsilon),\
                        prm.nonbond[i].sigma, math.fabs(epsilon14), sigma14, prm.nonbond[i].name) )
 
        

        self.write("Atoms\n")
        for i in range(0, psf.atomnum):
            nonbondid = nonbondtypes.index(psf.nonbond_parid[i])
            coor = crd.get_coor(i)
            self.write("%12d %8d %6d %7.3f %16.10f %16.10f %16.10f #%s"%(i+1, int(psf.atom[i].resid), nonbondid+1,\
                    psf.atom[i].charge, coor[0], coor[1], coor[2], psf.atom[i].atomstr ))
        
        if psf.bondnum > 0:
            self.write("\nBond Coeffs\n")
            bondtypes = list(bondtypes)
            bondtypes.sort()
            for i in range(0, len(bondtypes)):
                bondid = bondtypes[i]
                self.write("%12d %8.3f  %7.4f # %s"%(i+1, prm.bond[bondid].kb, prm.bond[bondid].b0," ".join(prm.bond[bondid].atom) ))

            self.write("\nBonds\n")
            for i in range(0, psf.bondnum):
                bondparid = bondtypes.index(psf.bond_parid[i]) + 1
                atomid1 = psf.bond[i][0]
                atomid2 = psf.bond[i][1]
                if atomid1 > atomid2:
                    atomid1 = psf.bond[i][1]
                    atomid2 = psf.bond[i][0]

                self.write("%12d %6d %9d %9d # %s %s"%(i+1, bondparid, atomid1, atomid2, \
                        psf.atom[atomid1-1].atomstr, psf.atom[atomid2-1].atomstr ))

        if psf.anglenum > 0:
            self.write("\nAngle Coeffs\n")
            angletypes = list(angletypes)
            angletypes.sort()
            constrained_angle = []
            for i in range(0, len(angletypes)):
                angleid = angletypes[i]
                if prm.angle[angleid].Kub < 0:
                    kub = 0
                    rub = 0
                else:
                    kub = prm.angle[angleid].Kub
                    rub = prm.angle[angleid].S0
                self.write("%12d %10.5f %10.5f %10.5f %10.5f # %s"%(i+1, prm.angle[angleid].Ktheta, \
                        prm.angle[angleid].Theta0, kub, rub, " ".join(prm.angle[angleid].atom)) )
                mass1 = prm.get_mass(prm.angle[angleid].atom[0])
                mass2 = prm.get_mass(prm.angle[angleid].atom[2])
                if int(mass1) == 1 and int(mass2) == 1:
                    constrained_angle.append(str(i+1))
            if len(constrained_angle) > 0:
                fpangle=open(lammpsangle, "w")
                fpangle.write("variable constraint_angletype  string \"%s\"\n"%(" ".join(constrained_angle)))
                fpangle.close()

            self.write("\nAngles\n")
            for i in range(0, psf.anglenum):
                angleparid = angletypes.index(psf.angle_parid[i]) + 1
                id1, id2, id3 = psf.angle[i]
                names = [psf.atom[j-1].atomstr for j in psf.angle[i]]
                names = " ".join(names)
                self.write("%12d %12d %12d %12d %12d # %s"%(i+1, angleparid, id1, id2, id3, names))

       
        if psf.dihenum > 0:
            self.write("\nDihedral Coeffs\n")
            for i in range(0, len(lammps_dihe_pars)):
                self.write("%12d %s"%(i+1, lammps_dihe_pars[i]))
            self.write("\nDihedrals\n")
            for i in range(0, len(lammps_dihe_lines)):
                self.write("%12d %s"%(i+1, lammps_dihe_lines[i]))

        if psf.imprnum > 0:
            self.write("\nImproper Coeffs\n")
            imprtypes = list(imprtypes)
            imprtypes.sort()
            for i in range(0, len(imprtypes)):
                imprid = imprtypes[i]
                self.write("%12d %8.3f %8.2f # %s"%(i+1, prm.impr[imprid].kpsi,\
                    prm.impr[imprid].psi0, " ".join(prm.impr[imprid].atoms) ) )
        
            self.write("\nImpropers\n")
            for i in range(0, psf.imprnum):
                imprparid = imprtypes.index(psf.impr_parid[i]) + 1
                atomnames = [psf.atom[j-1].atomname+"-"+psf.atom[j-1].atomtype for j in psf.impr[i]]
                atomnames = " ".join(atomnames)
                self.write("%12d %12d %12d %12d %12d %12d #%s"%(i+1, imprparid,\
                    psf.impr[i][0], psf.impr[i][1], psf.impr[i][2], psf.impr[i][3], atomnames))
        
        #CMAP
        if len(psf.cmap) > 0:
            self.write("\nCMAP\n")
            for i in range( 0, psf.cmapnum):
                type = [psf.atom[j-1].atomname + "-" + psf.atom[j-1].atomtype for j in psf.cmap[i]]
                type = " ".join(type)
                self.write("%12d %3d %12d %12d %12d %12d %12d #%s"%(i+1, psf.cmap_parid[i]+1, psf.cmap[i][0], psf.cmap[i][1], psf.cmap[i][2], psf.cmap[i][3], psf.cmap[i][7], type))

            if lammpscmap:
                fpcmap = open(lammpscmap, "w")
                for i in prm.cmap:
                    i.write(fpcmap)
                fpcmap.close()


def get_toppar(toppar):
    fp = open(toppar, "r")
    fs = []
    for line in fp.readlines():
        line.lstrip()
        line.strip("\n")
        if len(line) == 0:
            continue
        if line[0] in ["*", "!"]:
            continue
        linee = line.split()
        if len(linee) == 0:
            continue
        if linee[0].lower() == "stream":
            fs.append(linee[1])
            continue
        if linee[0].lower() == "open" and linee[1].lower() == "read":
            fs.append(linee[6])
    return fs

def main():
    des= """convert a charmm system to Desmond dms or LAMMPS data file"""
    parser = argparse.ArgumentParser(description=des, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-psf', action="store", dest="psf", required=True,\
            help="charmm psf file in xplor_ext format")
    parser.add_argument('-crd', action="store", dest="crd", required=False,\
            help="charmm crd file")
    parser.add_argument('-vel', action="store", dest="vel", required=False,\
            help="namd vel file to set the atom velocity; optional")
    parser.add_argument('-coor', action="store", dest="coor", required=False,\
            help="namd coor file to set the atom coordinates, overwrite -crd; optional")
    parser.add_argument('-t', action="store", dest="temperature", type=float, required=False,\
            help="simulation temperature in Kelvin, to generate velocity for each atom; skipped if -vel is used")
    parser.add_argument('-toppar', action="store", dest="toppar", required=True,\
            help="toppar.str from charmm-gui")
    parser.add_argument('-dms', action="store", dest="outputdms",  required=False,\
            help="output dms file")
    parser.add_argument('-lammps', action="store", dest="outputlammps",  required=False,\
            help="output lammps file")
    parser.add_argument('-lammpscmap', action="store", dest="outputlammpscmap",  required=False,\
            help="output CMAP file for lammps")
    parser.add_argument('-lammpsangle', action="store", dest="outputlammpsangle",  required=False,\
            help="restraints on angles for lammps")
    parser.add_argument('-boxx', action="store", dest="boxx", nargs=3, type=float, required=True,\
            help="box size x")
    parser.add_argument('-boxy', action="store", dest="boxy", nargs=3, type=float, required=True,\
            help="box size x")
    parser.add_argument('-boxz', action="store", dest="boxz", nargs=3, type=float, required=True,\
            help="box size x")
    parser.add_argument('-zcenter', action="store", dest="zcenter", type=float, required=True, default=0,\
            help="z center of the box")
    parser.add_argument('-posref', action="append", dest="posref", nargs=1, required=False,\
            help="crd file of the position restrained atoms (do not write unrestrained atoms to this file);\ncan be repeated for multiple files\nthe coordinates are the restraint target;\nfor equilibration in Desmond; optional")
    parser.add_argument('-posk', action="append", dest="posk", nargs=3, type=float, required=False,\
            help="restraint force constant fx, fy, and fz;\nfor equilibration in Desmond; optional")

    inputarg = parser.parse_args()
    
    psffile = inputarg.psf
    crdfile = inputarg.crd
    temperature = inputarg.temperature
    toppar = inputarg.toppar
    outputdms = inputarg.outputdms
    outputlammps = inputarg.outputlammps
    outputlammpscmap = inputarg.outputlammpscmap
    outputlammpsangle = inputarg.outputlammpsangle
    boxx = inputarg.boxx
    boxy = inputarg.boxy
    boxz = inputarg.boxz
    zcenter = inputarg.zcenter



    inputpsf = psf(psffile)
    if inputarg.crd:
        inputcrd = crd(crdfile)
        if inputpsf.atomnum != inputcrd.atomnum:
            raise Exception("number of atoms in psf and crd do not match: %d in psf %d in crd"%(inputpsf.atomnum, inputcrd.atomnum))    #namd coordinates

    if inputarg.coor:
        inputcrd = namdcoor(inputarg.coor)
        if inputcrd.natom != inputpsf.atomnum:
            raise Exception("atom number in %s and %s does not match"%(inputarg.coor, psffile))
 
    posref = inputarg.posref
    posk = inputarg.posk
    if posref and len(posref) != len(posk):
        raise Exception("number of positon restraint file and force constant set do not match")
    
    #prmfile = glob.glob(topparpath+"/*")
    #prmfile += glob.glob(topparpath+"/*.prm")
    prmfile = get_toppar(toppar)

    charmmprm = prm()
    for pfile in prmfile:
        if not os.path.isdir(pfile):
            charmmprm.read_prm(pfile)
    charmmprm.stat()
    inputpsf.get_parameter(charmmprm)
    
    boxsize = [boxx, boxy, boxz]

    #position restraints
    posrestraint = []
    if posref:
        for i in range(0, len(posref)):
            if os.path.isfile(posref[i][0]):
                cc = crd(posref[i][0])
                cc.add_posk(posk[i])
                posrestraint.append(cc)
            else:
                print posref[i][0], "does not exist, skip"

    #namd velocity
    velocity = []
    if inputarg.vel:
        vel = namdvel(inputarg.vel)
        velocity = vel.velocity
        if vel.natom != inputpsf.atomnum:
            raise Exception("atom number in %s and %s does not match"%(inputarg.vel, psffile))
    

    if outputdms != None:
        dmsfile = dms(outputdms)
        dmsfile.generate_dms(inputpsf, inputcrd, charmmprm, temperature, boxsize, posrestraint, velocity)
    if outputlammps != None:
        lammpsfile = lammps(outputlammps, "CHARMM-GUI lammps input")
        lammpsfile.generate_lammps(inputpsf, inputcrd, charmmprm, boxsize, outputlammpscmap, outputlammpsangle, zcenter)

if __name__ == "__main__":
    main()
