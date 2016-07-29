# -*- coding: utf-8 -*-
#pKcompare v2.3 by David - Shortened some stuff, made it faster by like 20%.
#Pulls the pK's from the with_cl and wo_cl pK.outs for a protein
#Returns the difference

from collections import OrderedDict,namedtuple
import sys
import numpy as np
import matplotlib.pyplot as plt

#Gets the names for all the residues
def nameget(prot):
    f = open("/Users/codyduell/Desktop/" + prot + "/with_cl/pK.out", "r")
    #f = open("C:/Users/judy/Documents/David/pKi.out","r")
    names = []
    for line in f:
        col = line.split()
        names.append(col[0])
    del names[0]
    f.close()
    return names

#Gets the pK values for with_cl
def linegetwith(prot):
    wicl = open("/Users/codyduell/Desktop/" + prot + "/with_cl/pK.out", "r")
    wipk = []
    for line in wicl:
        col = line.split()
        wipk.append(col[1])
    del wipk[0]
    #print wipk
    wicl.close()
    return wipk

#Gets the pK values for wo_cl
def linegetwithout(prot):
    wocl = open("/Users/codyduell/Desktop/" + prot + "/wo_cl/pK.out", "r")
    wopk = []
    for line in wocl:
        col = line.split()
        wopk.append(col[1])
    del wopk[0]
    #print wopk
    wocl.close()
    return wopk

#Changes pK values for with_cl, removes "<" and ">" signs.
def normalwith(prot):
    wipk = linegetwith(prot)
    for index in range(len(wipk)):
        if wipk[index] == ">14.0":
            wipk[index] = "14.0"
        elif wipk[index] == "<0.0":
            wipk[index] = "0.0"
    #print wipk
    return wipk
    
#Changes pK values for wo_cl, removes "<" and ">" signs.
def normalwithout(prot):
    wopk = linegetwithout(prot)
    for index in range(len(wopk)):
        if wopk[index] == ">14.0":
            wopk[index] = "14.0"
        elif wopk[index] == "<0.0":
            wopk[index] = "0.0"
    #print wopk
    return wopk

#Calls on each of the previous functions to create a dictionary of values, arranged by (Residue Name, pK difference)
#Prints said dictionary and returns residues with a pK difference of 4(?) or more.
#Note: the difference is generated as pKwith - pKwithout.
def compare(prot):
    wipk = normalwith(prot)
    wopk = normalwithout(prot)
    difference = []
    for line in range(len(wipk)):
        difference.append(round(float(wipk[line])-float(wopk[line]),3))
    #print difference
    names = nameget(prot)
    if len(names) != len(difference):
        print "Well, something went wrong. The number of names isn't equal to the number of results. Consult the PDB or the script for errors."
        sys.exit()
    diffs = OrderedDict(zip(names,difference))
    return diffs
    for item in range(len(diffs)):
        print diffs.items()[item]
    print "Notable differences:"
    for item in range(len(diffs)):
        if abs(diffs.items()[item][1]) >= 4:
            print diffs.items()[item] 

#Gets all of the names of the residues from the pK.outs
def reslist(prot):
    diffs = compare(prot)
    reslist = diffs.keys()
    return reslist

#Gets the coordinates for the chlorides
def getcl(prot):
    pdb = open("/Users/codyduell/Desktop/" + prot + "/prot.pdb", "r")
    cl = []
    for line in pdb:
        if "HETATM" in line and "CL" in line:
            tmp = []
            col = line.split()
            if len(col[4]) == 5:
                col.insert(4, col[4][0:1])
                col.insert(5, col[5][1:])
                col.pop(6)
            tmp.extend(list([float(col[6]),float(col[7]),float(col[8])]))
            cl.append(tmp)
    pdb.close()
    return cl   

#Gets the lines for all the residues and chlorides
def getlines(prot):
    pdb = open("/Users/codyduell/Desktop/" + prot + "/prot.pdb", "r")
    atoms = []
    for line in pdb:
        col = line.split()
        if col[0] == "ATOM":
            atoms.append(line)
        if "HETATM" in line and "CL" in line:
            atoms.append(line)
    for line in range(len(atoms)):
        col = atoms[line].split()
        if len(col[2]) >= 4:
            col.insert(2,col[2][0:4])
            col.insert(3,col[3][4:])
            col.pop(4)
        if len(col[4]) == 5:
            col.insert(4, col[4][0:1])
            col.insert(5, col[5][1:])
            col.pop(6)
        if len(col) != 12:
            print col
            print "I'm just a computer, but I'm pretty sure this is wrong."
            sys.exit()
        atoms[line] = col
    pdb.close()
    return atoms

#Gets all the coordinates
def matchresidues(prot):
    resl = reslist(prot)
    atoms = getlines(prot)
    rcoord = []
    for col in atoms:
        for res in resl:
            xyz = []
            if res[0:3] == col[3] and col[2] == "CA" and int(res[5:9]) == int(col[5]):
                resl.remove(res)
                xyz.extend([float(col[6]),float(col[7]),float(col[8])])
                rcoord.append(xyz)
            if res[1:3] == col[3] and col[0] == "HETATM" and int(res[5:9]) == int(col[5]):
                resl.remove(res)
                xyz.extend([float(col[6]),float(col[7]),float(col[8])])
                rcoord.append(xyz)
    res = namedtuple('Result',['first','second'])
    r = res(rcoord,resl)
    return r

#Figures out which residue is the n-term
def nter(prot):
    nt = matchresidues(prot)[1][0]
    atoms = getlines(prot)
    for atom in atoms:
        if (int(nt[5:9]) == int(atom[5])) and atom[2] == "CA":
            nterm = str(nt) + "= " + atom[3] + "_" + atom[5]
            break
    return nterm

#Figures out which residue is the c-term
def cter(prot):
    ct = matchresidues(prot)[1][1]
    atoms = getlines(prot)
    for atom in atoms:
        if (int(ct[5:9]) == int(atom[5])) and atom[2] == "CA":
            cterm = str(ct) + "= " + atom[3] + "_" + atom[5]
            break
    return cterm

#Figures out all of the distances, returns the shortest one
def measure(prot):
    dist = []
    rcoord = matchresidues(prot)[0]
    ccoord = getcl(prot)
    for m in range(len(ccoord)):
        c = np.asarray(ccoord[m])
        tempdist = []
        for n in range(len(rcoord)):
            r = np.asarray(rcoord[n])
            tempdist.append(float(np.linalg.norm(r-c)))
        dist.append(tempdist)
    for i in range(len(dist)):
        for k in range(len(dist[i])):
            dist[i][k] = round(dist[i][k],2)
    dist = zip(*dist)
    shortest = []
    for i in dist:
        shortest.append(min(list(i)))
    return shortest
    #line 203-208 can be a separate function
    
#Combines it all and prints it out
def compiled():
    prot = raw_input("Prot name:")
    diffs = compare(prot)
    shortest = measure(prot)
    nt = matchresidues(prot)[1][0]
    ct = matchresidues(prot)[1][1]
    nname = nter(prot)
    cname = cter(prot)
    threes = []
    for item in range(len(diffs.items())):
        if not (diffs.keys()[item] == nt or diffs.keys()[item] == ct):
            threes.append(list(diffs.items()[item]))
    if len(threes) != len(shortest):
        print "NOPE"
        sys.exit()
    for distance in range(len(shortest)):
        threes[distance].append(shortest[distance])
    print nname
    print cname
    print "#, Res, difference, dist from nearest cl"
    for i in range(len(threes)):
        print "{0:3d}".format(i+1), threes[i]
    print "\n", "Notable differences:" 
    for item in range(len(threes)):
        if abs(threes[item][1]) >= 4:
            print threes[item]
    return threes
    
def plot():
    threes = compiled()
    diffs = []
    dists = []
    for i in range(len(threes)):
        diffs.append(abs(threes[i][1]))
        dists.append(abs(threes[i][2]))
    #print diffs
    #print dists
    plt.scatter(diffs,dists)
    plt.gca().set_xlim(left=-0.27)
    plt.gca().set_ylim(bottom=-0.27)
    plt.show()
plot()

