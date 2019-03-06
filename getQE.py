#!/usr/bin/env pythonw
'''
Program header:
    Mingfei Zhang
run with "pythonw" as framework
code to verify quasi-entropy on binary alloy and complex structure binary alloy
python3 code
need to install ase first
https://wiki.fysik.dtu.dk/ase/
'''
import os
import sys
import numpy as np
from ase.io import read
from ase import neighborlist
import matplotlib.pyplot as plt
from scipy import spatial

def get_cos_dist(dataSetI, dataSetII):
    return spatial.distance.cosine(dataSetI, dataSetII)

class QEntropy(object):

    def get_atom_F_A_B(self, atm, types, nBins=200):
        rdf_Atom_all=np.zeros((len(types),nBins))
        #print(rdf_Atom_all.shape)
        #get index of element types of interest
        k=0
        for iType in types:
            typeIndex=0
            for j in range(len(self.listTypeName)):
                if iType == self.listTypeName[j]:
                    typeIndex = j
            typeBStart = 0
            typeBEnd = 0
            for i in range(len(self.listTypeNum)):
                if i < typeIndex:
                    typeBStart += self.listTypeNum[i]
            typeBEnd = typeBStart + self.listTypeNum[typeIndex]
    
            atoms_A = self.atoms[atm:atm+1]
            atoms_B = self.atoms[typeBStart:typeBEnd]
            atoms_new = atoms_A + atoms_B
            d = neighborlist.neighbor_list('d', atoms_new, self.rMax)
            dr = self.rMax/nBins
            edges = np.arange(0., self.rMax + 1.1 *dr, dr)
            h, binEdges = np.histogram(d, edges)
            rho = len(atoms_new) / self.atoms.get_volume() 
            factor = 4./3. * np.pi * rho * len(atoms_new)
            rdf = h / (factor * (binEdges[1:]**3 - binEdges[:-1]**3)) 

            rdf_Atom_all[k,:] = rdf[1:]-1 #skip first one in case 0 dist
            k+=1
        #print(rdf_Atom_all)
        #plt.plot(binEdges[1:nBins+1], rdf_Atom_all[0])
        #plt.plot(binEdges[1:nBins+1], rdf_Atom_all[1])
        #plt.savefig('RDF/rdf_'+str(atm)+'.pdf')
        #plt.close()
        return(rdf_Atom_all)

    def cal_QE(self):
        S_str = 0.0
        l = 0
        for iType in self.listTypeName:
            typeIndex=0
            for j in range(len(self.listTypeName)):
                if iType == self.listTypeName[j]:
                    typeIndex = j
            typeAStart = 0
            typeAEnd = 0
            for k in range(len(self.listTypeNum)):
                if k < typeIndex:
                    typeAStart += self.listTypeNum[k]
            typeAEnd = typeAStart + self.listTypeNum[typeIndex]
            tmp = 0.0
            count = 0
            for i in range(typeAStart, typeAEnd-1):
                for j in range(i+1, typeAEnd):
                    tmp += np.log(1.0 - get_cos_dist(self.F[i][l],self.F[j][l]))
                    count+=1
            tmp/=count
            tmp = tmp *self.listTypeNum[l]/self.tot_num
            S_str += tmp
        l+=1
        S_str = -S_str
        return(S_str)

    def get_All_F(self):
        for i in range(self.tot_num):
            self.F[i,:,:] = self.get_atom_F_A_B(i, self.listTypeName, self.nBins)
        return

    def __init__(self, inFile, rMax = 6.0, nBins = 100):
        self.inFile = inFile
        print("working on: %s" % self.inFile)
        self.rMax = rMax
        self.nBins = nBins
        with open(inFile,'r') as fin:
            lines=fin.readlines()
            self.listTypeName=[str(s) for s in lines[5].split()]
            self.listTypeNum=[int(s) for s in lines[6].split()]
        self.tot_num=int(sum(self.listTypeNum))
        self.atoms = read(self.inFile)
        self.F = np.zeros((self.tot_num, len(self.listTypeName), self.nBins))
        self.get_All_F()
        print(self.cal_QE())
        return
    
    
if __name__ == "__main__":
    if (len(sys.argv) != 2):
        print("usage:\ngetQE.py <vasp5 formate file name>\n")
    else:
        inFile = sys.argv[1]
        cnf = QEntropy(inFile)
