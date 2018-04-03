# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 16:19:23 2018

@author: HuangMing
"""

import os
import pandas as pd
import itertools
import numpy as np

import AtomClass as Atom
import MoleculeClass as Mol
import SystemClass as SysClass

def splitString(s):
    string = []
    digit = []
    for i in s:
        if i.isdigit():
            digit.append(i)
        else:
            string.append(i)
    string = ''.join(string)
    dig = ''.join(digit)
    
    return dig, string

def InputSysInfo(name, num, size):
    system = SysClass.SystemInfo()
    system.setName(name)
    system.setNum(num)
    system.setSize(size)

def AtomInfoUpdate(atom, idx, Index=False):
    if (Index):
        atom.setglobalIndex(idx)
                    
def AtomInfoInput(atom, line):
    str1 = splitString(line[0])
    molNum = str1[0]
    molName = str1[1]
    atomName = line[1]
    atomIndex = line[2]
    atomPos = [round(float(line[3]), 3), round(float(line[4]), 3), round(float(line[5]), 3)]
    atom.setmolNum(molNum)   
    atom.setmolName(molName)
    atom.setatomName(atomName)
    atom.setglobalIndex(atomIndex)
    atom.setPos(atomPos)
    return atom

def DelAtoms(atomName, mol): #Doesn't work, needs to be fixed, suppose should because the mol doesn't change the mollist content
    print('del atoms: ', atomName)
    for atom in mol.getAtoms():
        if atom.getatomName == atomName:
            index = atom.getlocalIndex()
            mol.pop(index)
    
def MolInfoInput(index, name, atomList):
    mol = Mol.MoleculeInfo()
    mol.setIndex(index)
    mol.setName(name)
    for i in range(len(atomList)):
        mol.setAtoms(atomList[i])
        
    return mol

def MolInfoUpdate(mol, resName, Name=False):
    if Name:
        mol.setName(resName)
        atoms = mol.getAtoms()
        for atom in atoms:
            atom.setmolName(resName)
            name = atom.getatomName()
            if name == 'C1':
                atom.setatomName('C')
            elif name == 'N1':
                atom.setatomName('N')
                
def GetTopMol(atomList):
    index = atomList[0].getmolNum()
    name = atomList[0].getmolName()
    i = 0
    localIndex = 0
    while (i < len(atomList)):
        idx = atomList[i].getmolNum()
        atomList[i].setlocalIndex(localIndex)
        if idx != index:
            mol = MolInfoInput(index, name, atomList[0:i])
            return i, mol
        i += 1
        localIndex += 1
    index = atomList[i-1].getmolNum()
    name = atomList[i-1].getmolName()
    mol = MolInfoInput(index, name, atomList[0:i])
    return i, mol

def SplitAtom(baseList):
    atomsList = []
    for i in range (len(baseList)):
        atom = Atom.AtomsInfo()
        atom = atomsList.append(AtomInfoInput(atom, baseList.iloc[i]))
    
    return atomsList

def Atom2Mol(atomsList):
    molList = []
    i = 0
    num = 1
    while i < len(atomsList):
        info = GetTopMol(atomsList[i:])
        mol = info[1]
        molList.append(mol)
        num += 1
        i += info[0]
    return molList

def CheckReactAtoms(atomName, atomsList):
    for i in range (len(atomsList)):
        name = atomsList[i].getatomName()
        if name == atomName:
            return True
    return False

def GetReactAtoms(atomName, atomsList):
    index = []
    for i in range (len(atomsList)):
        name = atomsList[i].getatomName()
        if name == atomName:
#            idx = atomsList[i].getglobalIndex()
            index.append(i)
    return index

def SplitMonCro(monR, croR, molList):
    monReactiveAtoms = []
    croReactiveAtoms = []
    for i in range (len(molList)):
        atoms = molList[i].getAtoms()
        if (CheckReactAtoms(monR, atoms)):           
            idxMon = GetReactAtoms(monR, atoms)
            monReactiveAtoms.append([i, idxMon])
        elif (CheckReactAtoms(croR, atoms)):
            idxCro = GetReactAtoms(croR, atoms)
            croReactiveAtoms.append([i, idxCro])
        else:
            pass
#            print('No specified atoms in the list, please check your input!')
    return monReactiveAtoms, croReactiveAtoms

def CalDist(aPos, bPos):
    x2 = np.power(float(aPos[0])-float(bPos[0]),2)
    y2 = np.power(float(aPos[1])-float(bPos[1]),2)
    z2 = np.power(float(aPos[2])-float(bPos[2]),2)
    dist = np.sqrt(x2+y2+z2)
    return dist

def GetMonCroList(monR, croR, molList):
    reactLists = SplitMonCro(monR, croR, molList)
    monAtoms = []
    croAtoms = []
    monList = reactLists[0]
    croList = reactLists[1]
    
    for mon in monList:
        atom = GetAtom(mon, molList)
        monAtoms.append(atom)

    for cro in croList:
        atom = GetAtom(cro, molList)
        croAtoms.append(atom)
    list(itertools.chain.from_iterable(monAtoms))
    list(itertools.chain.from_iterable(croAtoms))
    
    return monAtoms, croAtoms

def Crosslink(monR, croR, molList, cutoff):
    initList = GetMonCroList(monR, croR, molList)
    croAtoms = initList[1]
    outputName = 'tst.top'
#    os.remove('tst.top')
    for croAtom in croAtoms:
        initList = GetMonCroList(monR, croR, molList)
        monAtoms = initList[0]
        posCRO = croAtom[0].getPos()

        for monAtom in monAtoms:
            posMON = monAtom[0].getPos()
            Genbond(monAtom[0], croAtom[0], molList, cutoff, outputName) 
    return molList

def Genbond(atom1, atom2, molList, cutoff, outputName):
    reactedMON = 'MOO'
    reactedCRO = 'COO'
    atom1_Idx = atom1.getglobalIndex()
    atom1_localIdx = atom1.getlocalIndex()
    atom2_Idx = atom2.getglobalIndex()
    atom2_localIdx = atom2.getlocalIndex()
    
    mol1_Idx = int(atom1.getmolNum())
    mol2_Idx = int(atom2.getmolNum())
    dist = CalDist(atom1.getPos(), atom2.getPos())
    cond1 = Criteria1(dist, cutoff)
#    print(mol2_Idx)
#    print(len(molList))
    if cond1:
        mol1 = molList[mol1_Idx - 1]
        mol2 = molList[mol2_Idx - 1]
        MolInfoUpdate(mol1, reactedMON, Name=True)
        MolInfoUpdate(mol2, reactedCRO, Name=True)
        DelAtoms('H1',mol1)
        f = open(outputName, 'a')
        str1 = "Idx1: {:>6},    Idx2: {}\n".format(atom1_Idx, atom2_Idx)
        f.write(str1)
        f.close()
#    f = open('tst.top', 'a')
#    
#    str1 = "Idx: {}, localIdx: {}, molIdx: {}\n".format(atom1_Idx, atom1_localIdx, mol1_Idx)
#    f.write(str1)
#    f.close()
    
#    return monList, croList
#    AtomInfoUpdate(molList[mol1_Idx], )

def Criteria1(dist, cutoff): 
    if dist < cutoff:
        return True
    else:
        return False

def GetAtom(Index, molList):
    atoms = []
    molIndex = Index[0]
    atomIndex = Index[1]
    for idx in atomIndex:    
        atom = molList[molIndex].getAtoms()[idx]
        atoms.append(atom)
    return atoms

def ExportGRO(info, outputName, molList):
    index = 1
    for i in range (len(molList)):
        length = len(molList[i].getAtoms())
        atoms = molList[i].getAtoms()
        for ii in range (len(atoms)):
            idx = index + ii
            atom = atoms[ii]
            AtomInfoUpdate(atom, idx, Index=True)
        index += length
#    for atom in atomsList:
#        str1 = atom.outputInfo()
#        print(str1)
    f = open(outputName, 'w')
    molName = info[0] + '\n'
    molNum = info[1] + '\n'
    molSize = info[2] + '\n'
    
    f.write(molName)
    f.write(molNum)
    for mol in molList:
        atoms = mol.getAtoms()
        for atom in atoms:
            str1 = atom.outputInfo() + '\n'
#            print(str1)
            f.write(str1)
    f.write(molSize)
    f.close()
    
filename = 'min.gro'
outputName = 'tst.gro'
monR = 'C1'
croR = 'N1'
cutoff = 1.

#df = pd.read_json(filename, sep='\n', header=None)
df = pd.read_csv(filename, sep='\n', header=None)
baseList = df.iloc[2:-1].reset_index(drop=True)[0].str.split()
molName = df.iloc[0][0]
molNum = df.iloc[1][0]
molSize = df.iloc[-1][0]
info = [molName, molNum, molSize]

atomsList = SplitAtom(baseList)
molList = Atom2Mol(atomsList)
#ReactAtoms = SplitMonCro(monR, croR, molList)
a = Crosslink(monR, croR, molList, cutoff)
ExportGRO(info, outputName, a)
#a = df.iloc[2][0]
#b = splitString(a)