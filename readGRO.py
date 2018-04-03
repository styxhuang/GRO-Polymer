# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 16:19:23 2018

@author: HuangMing
"""

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

def DelAtoms(index, molList):
    molList.pop(index)
    
def MolInfoInput(index, name, atomList):
    mol = Mol.MoleculeInfo()
    mol.setIndex(index)
    mol.setName(name)
    for i in range(len(atomList)):
        mol.setAtoms(atomList[i])
        
    return mol

def MolInfoUpdate():
    pass

def GetTopMol(atomList):
    index = atomList[0].getmolNum()
    name = atomList[0].getmolName()
    i = 0
    while (i < len(atomList)):
        idx = atomList[i].getmolNum()
        if idx != index:
            mol = MolInfoInput(index, name, atomList[0:i])
            return i, mol
        i += 1
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
#        print('tst-2:', i)
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
            print('No specified atoms in the list, please check your input!')
    return monReactiveAtoms, croReactiveAtoms

def CalDist(aPos, bPos):
    x2 = np.power(float(aPos[0])-float(bPos[0]),2)
    y2 = np.power(float(aPos[1])-float(bPos[1]),2)
    z2 = np.power(float(aPos[2])-float(bPos[2]),2)
    dist = np.sqrt(x2+y2+z2)
    return dist

def Crosslink(reactLists, molList, cutoff):
    monAtoms = []
    croAtoms = []
    monList = reactLists[0]
    croList = reactLists[1]
#    print(monList)
    for mon in monList:
        atom = GetAtom(mon, molList)
        monAtoms.append(atom)

        
    for cro in croList:
        atom = GetAtom(cro, molList)
        croAtoms.append(atom)
    
    list(itertools.chain.from_iterable(monAtoms))
    list(itertools.chain.from_iterable(croAtoms))
    
    for croAtom in croAtoms:
        posCRO = croAtom[0].getPos()
        for monAtom in monAtoms:
            posMON = monAtom[0].getPos()
            dist = CalDist(posCRO, posMON)
            if dist < cutoff:
                print(monAtom[0].getglobalIndex())
                print("bond!")
    return monAtoms, croAtoms
    
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
cutoff = 3.

#df = pd.read_json(filename, sep='\n', header=None)
df = pd.read_csv(filename, sep='\n', header=None)
baseList = df.iloc[2:-1].reset_index(drop=True)[0].str.split()
molName = df.iloc[0][0]
molNum = df.iloc[1][0]
molSize = df.iloc[-1][0]
info = [molName, molNum, molSize]

atomsList = SplitAtom(baseList)
molList = Atom2Mol(atomsList)
ReactAtoms = SplitMonCro(monR, croR, molList)
a = Crosslink(ReactAtoms, molList, cutoff)
ExportGRO(info, outputName, molList)
#a = df.iloc[2][0]
#b = splitString(a)