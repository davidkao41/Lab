#!/usr/bin/env/python
#v9 Fixed a problem with results being displayed twice. Fixed a small error with cutoffs.
# Minor fix for more reliable line-pulling in atomherding().
#Added comments for future updates
# %timeit (command line) to time it | >>>%timeit ("filepath")
#Can: compare 2 residues, compare CL with another residue
#use full if you want to just see all distances between the two atoms between all residues
from collections import OrderedDict
from math import sqrt
import sys
atoms = []
findres1 = []
findres2 = []
everything = []
first= []
second = []
chl = []
temp = []
mytwo = []
#Crossroads asks you for the variables you want. Findatom looks for the atom,
#and the res's bring up the certain residues. Full asks if you want every distance
#between every residue. Cutoff returns the distances under the cutoff.
def crossroads():
    global findatom1,findatom2,full, firstres, secondres, cutoff,mytwo
    findatom1 = sys.argv[1]#caps
    findatom2 = sys.argv[2]#caps
    full = sys.argv[3]#type y/n
    if full != "y":
        firstres = raw_input("3 letter res code ->") #use 3 CAPITAL letters, like GLN
        secondres = raw_input("3 letter res code -->")#if you pick CL, this will be ignored
        cutoff = raw_input("cutoff --->")#any number is ok
        mytwo.append(firstres)
        mytwo.append(secondres)
#Pulls the information for each atom from the pdb.
def atomherding():
    with open("C:/Users/Kao_2/Downloads/lab/2g55.pdb") as pdb:
        for line in pdb:
            col = line.split()
            if col[0] == "ATOM":
                atoms.append(line)
#Pulls the information for each CL from the pdb.
def atomherdingcl():
    with open("C:/Users/Kao_2/Downloads/lab/2g55.pdb") as pdb:
        for line in pdb:
            if "HETATM" in line:
                temp.append(line)
        new = list(OrderedDict.fromkeys(temp))
        for line2 in new:
            colhetatm = line2.split()
            if colhetatm[2] == "CL":
                chl.append(colhetatm)
#Checks to see if the lines are correctly formatted -> checks to see if each
#line has 12 columns. Includes a linefixer for a common formatting error in
#column 3 and 4.
def pdbisright():
    for a in range(len(atoms)):
        column = atoms[a].split()
        if len(column[2]) >= 4:
            column.insert(2,column[2][0:4])
            column.insert(3,column[3][4:])
            column.pop(4)
            print column,"\n"
        if len(column) != 12:
            print "PDB is wrong"
            print column
            sys.exit()
#Makes sure the CL line is correctly formatted.
def pdbisrightcl():
    for b in range(len(chl)):
        if len(chl[b]) != 12:
            print "chl is wrong"
            print chl[b]
            sys.exit()
#puts all of the atoms with the element names you want, and pulls them into a list
def atomculling():
    for atom in range(len(atoms)):
        column = atoms[atom].split()
        if column[2] == findatom1:
            findres1.append(atoms[atom])
        if column[2] == findatom2:
            findres2.append(atoms[atom])
#puts all of the atoms with the element name you want, and pulls them into a list
def atomcullingcl():
    for atom in range(len(atoms)):
        column = atoms[atom].split()
        if column[2] == findatom2:
            findres2.append(atoms[atom])
#calculates your distance between the two atoms
def distance():
    for atom1 in range(len(findres1)):
        column1 = findres1[atom1].split()
        for atom2 in range(len(findres2)):
            column2 = findres2[atom2].split()
            deltax = abs(float(column1[6])-float(column2[6]))
            deltay = abs(float(column1[7])-float(column2[7]))
            deltaz = abs(float(column1[8])-float(column2[8]))
            dist = sqrt((deltax**2) + (deltay**2) + (deltaz**2))
            Results = []
            Results.append([int(column1[1]),column1[3],column1[2]])
            Results.append([int(column2[1]),column2[3],column2[2]])
            Results.append(dist)
            everything.append(Results)
#calculates the distance between the cl and the other atom.
def distancecl():
    for atom1 in range(len(chl)):
        for atom2 in range(len(findres2)):
            column2 = findres2[atom2].split()
            deltax = abs(float(chl[atom1][6])-float(column2[6]))
            #Add a preliminary cutoff. Advice needed?
            deltay = abs(float(chl[atom1][7])-float(column2[7]))
            deltaz = abs(float(chl[atom1][8])-float(column2[8]))
            #Remove the sqrt
            dist = sqrt((deltax**2) + (deltay**2) + (deltaz**2))
            Results = []
            Results.append([int(chl[atom1][1]),chl[atom1][3],chl[atom1][2]])
            Results.append([int(column2[1]),column2[3],column2[2]])
            Results.append(dist)
            everything.append(Results)
#For printing all distances between the atoms without regard for residues
def fullprint():
    if full == "y":
        for b in range(len(everything)):
            print everything[b]
        sys.exit()
#Sorts results and prints out the ones with the two residues you entered and the cutoff.
def residuecom():
    for h in range(len(everything)):
        for i in range(2):
            if firstres in everything[h][i]:
                first.append(everything[h])
            if secondres in everything[h][i]:
                first.append(everything[h])
    for j in range(len(first)):
        if (first[j][0][1] in mytwo) and (first[j][1][1] in mytwo):
            if first[j] not in second:
                second.append(first[j])
    for q in range(len(second)):
        if second[q][2] < float(cutoff):
            print second[q]
#Sorts results and prints out the ones with the residue you entered and the cutoff.
def residuecomcl():
    for w in range(len(everything)):
        for y in range(2):
            if firstres in everything[w][y]:
                first.append(everything[w])
    for q in range(len(first)):
        if first[q][2] < float(cutoff):
            print first[q]
#Actual execution
def execl():
    atomherding()
    atomherdingcl()
    pdbisright()
    pdbisrightcl()
    atomcullingcl()
    distancecl()
    fullprint()
    residuecomcl()
def exen():
    atomherding()
    pdbisright()
    atomculling()
    distance()
    fullprint()
    residuecom()
def doit():
    crossroads()
    if findatom1 == "CL":
        execl()
    else:
        exen()
doit()
