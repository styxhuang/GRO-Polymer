# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 16:19:23 2018

@author: HuangMing
"""

import os
import pandas as pd
import itertools
import numpy as np
import random

import AtomClass as Atom
import MoleculeClass as Mol
import SystemClass as SysClass

dictBond = {
        'C-N': '6, 0.14700, 268278.'
        }

dictAngle = { #TODO: keep update angle coefficient
        'N-C-N': '1, 110.380, 553.9616',
        'C-O-C': '1, 117.600, 522.1632',
        'C-C-N': '1, 110.380, 553.9616',
        'C-N-C': '1, 110.900, 535.5520',
        'N-C-C': '1, 110.380, 553.9616',
        'H-C-N': '1, 109.920, 413.3792',
        'C-N-H': '1, 109.920, 394.1328'
        }

dictDihedral = { #Has little different on the 3 column, some combination is 2 and 3'
        'C-C-O-C': '1, 180.000, 3.766, 2',
        'C-O-C-C': '1, 180.000, 4.602, 2',
        'C-N-C-C': '1, 180.000, 2.008, 2',
        'C-C-N-C': '1, 180.000, 2.008, 2'
        }

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

def DelAtoms(atomName, mol):
    for atom in mol.getAtoms():
        if atom.getatomName() == atomName:
            index = atom.getlocalIndex()
            mol.getAtoms().pop(index)
#    return mol

def DelHydrogen(idx, mol):
    mol.getAtoms().pop(idx+1)
    atoms = mol.getAtoms()
    for i in range (len(atoms)):
        atoms[i].setlocalIndex(i)
        
#    for atom in mol.getAtoms():
#        index = atom.getlocalIndex()
#        print('localIdx: ', index)

def MolInfoInput(index, name, atomList):
    mol = Mol.MoleculeInfo()
    mol.setIndex(str(int(index) + 1)) #this mol index start from 0
    mol.setName(name)
    for i in range(len(atomList)):
        mol.setAtoms(atomList[i])
        
    return mol

def MolInfoUpdate(mol, resName, localIdx, atomBR, atomAR, Name=False):
    if Name:
        if localIdx != 0: #Assume if not the 1st atom reacted, then the residue name will add 'R' at the end
            resName = resName + 'R'
        mol.setName(resName)
        atoms = mol.getAtoms()
        print('localIdx', localIdx)
        atoms[localIdx].setatomName(atomAR)
        for atom in atoms:
            atom.setmolName(resName)        
                
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
            index.append(i)
    return index

def SplitMonCro(monR, croR, molList):
    monReactiveAtoms = []
    croReactiveAtoms = []
    print('molList', len(molList))
    for i in range (len(molList)):
        atoms = molList[i].getAtoms()
        if (CheckReactAtoms(monR, atoms)):           
            idxMon = GetReactAtoms(monR, atoms)
#            print('idxMon', idxMon)
            for idx in idxMon:
                monReactiveAtoms.append([i, idx])
        elif (CheckReactAtoms(croR, atoms)):
            idxCro = GetReactAtoms(croR, atoms)
#            print('idxCro', idxCro)
            for idx in idxCro:
                croReactiveAtoms.append([i, idx])
        else:
            pass
            print('No specified atoms in the list, please check your input!')
    print('monReactiveAtoms', len(monReactiveAtoms))
    print('croReactiveAtoms', len(croReactiveAtoms))
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
    print('monList', monList)
    croList = reactLists[1]
    print('croList', croList)
    for mon in monList:
        atom = GetAtom(mon, molList)
        monAtoms.append(atom)

    for cro in croList:
        atom = GetAtom(cro, molList)
        croAtoms.append(atom)
        
#    print('monAtoms', monAtoms)
#    print('croAtoms', croAtoms)
#    list(itertools.chain.from_iterable(monAtoms))
#    list(itertools.chain.from_iterable(croAtoms))
    
    return monAtoms, croAtoms

def Crosslink(monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, molList, cutoff, bondCycle):
    pairs = []
    initList = GetMonCroList(atom1BR, atom2BR, molList)
    croAtoms = initList[1]
    outputName = 'tst.top'
    if os.path.isfile(outputName):
        os.remove(outputName)
    bond = 0
    random.shuffle(croAtoms)
    for croAtom in croAtoms:
        initList = GetMonCroList(atom1BR, atom2BR, molList) #Since mon list needs to update each loop
#        print('initList: ', initList)
        monAtoms = initList[0]
        random.shuffle(monAtoms)

        for monAtom in monAtoms:
            pair = Genbond(monAtom, croAtom, molList, cutoff, outputName, monAR, croAR, atom1BR, atom1AR, atom2BR, atom2AR)  
            print('pair: ', pair)
            if pair != 0:
                pairs.append(pair)
                break
        bond += 1
        if bond >= bondCycle:
            break

    return pairs, molList

def Genbond(atom1, atom2, molList, cutoff, outputName, monAR, croAR, atom1BR, atom1AR, atom2BR, atom2AR):
    reactedMON = monAR
    reactedCRO = croAR

    mol1_Idx = int(atom1.getmolNum())
    mol2_Idx = int(atom2.getmolNum())
    local_Idx1 = atom1.getlocalIndex()
    local_Idx2 = atom2.getlocalIndex()
    
    dist = CalDist(atom1.getPos(), atom2.getPos())
    cond1 = Criteria1(dist, cutoff)
    
    if cond1:      
        mol1 = molList[mol1_Idx - 1]
        mol2 = molList[mol2_Idx - 1]
        MolInfoUpdate(mol1, reactedMON, local_Idx1, atom1BR, atom1AR, Name=True)
        MolInfoUpdate(mol2, reactedCRO, local_Idx2, atom2BR, atom2AR, Name=True)
        
        idx_1 = atom1.getlocalIndex()
        idx_2 = atom2.getlocalIndex()
        
        DelHydrogen(idx_1, mol1)
        DelHydrogen(idx_2, mol2)
        
        pair = [[mol1_Idx, atom1.getlocalIndex()],[mol2_Idx, atom2.getlocalIndex()]]
#        print('atom1_Idx_del: ', atom1.getglobalIndex())
#        f = open(outputName, 'a')
#        str1 = "{:>8}{:>6}{:>4}{:>12}{:>12}\n".format(atom1_Idx, atom2_Idx, 6, 0.14700, 268278.)
#        f.write(str1)
#        f.close()
        return pair
    return 0

def Criteria1(dist, cutoff): 
    if dist < cutoff:
        return True
    else:
        return False

def Criterial2():
    pass
def GetAtom(Index, molList):
    molIndex = Index[0]
    atomIndex = Index[1]
    
    atom = molList[molIndex].getAtoms()[atomIndex]
    return atom

def GetInfo(info, category='bond'):
    if category == 'bond':
        if info in dictBond:
            coeff = dictBond[info]
            return coeff
    else:
        print("Don't have this info")
    
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
    f = open(outputName, 'w')
    molName = info[0] + '\n'
    ##############################
    num = 0
    for mol in molList:
        atoms = mol.getAtoms()
        for atom in atoms:
            num += 1
    molNum = str(num) + '\n'
    ##############################
    molSize = info[2] + '\n'
    
    f.write(molName)
    f.write(molNum)

    for mol in molList:
        atoms = mol.getAtoms()
        for atom in atoms:
            str1 = atom.outputInfo() + '\n'
            f.write(str1)
    f.write(molSize)
    f.close()

def HasNumbers(inString):
    return any(char.isdigit() for char in inString)

def ExtraLetterCheck(inString):
    string = ['C', 'N', 'O']
    return [char for char in inString if char in string]

def GetBondInfo(molInfo, molList):
    atom1Info = molInfo[0]
    atom2Info = molInfo[1]
    
    mol1_idx = atom1Info[0] - 1
    mol2_idx = atom2Info[0] - 1

    atom1 = molList[mol1_idx].getAtoms()[atom1Info[1]]
    atom2 = molList[mol2_idx].getAtoms()[atom2Info[1]] 

    idx_1 = atom1.getglobalIndex()
    idx_2 = atom2.getglobalIndex()
    
    name1 = atom1.getatomName()
    name2 = atom2.getatomName()
        
#    if HasNumbers(name1):
#        name1 = ''.join([i for i in name1 if not i.isdigit()])
#    if HasNumbers(name2):
#        name2 = ''.join([i for i in name2 if not i.isdigit()])   
    name1 = ''.join(ExtraLetterCheck(name1))
    name2 = ''.join(ExtraLetterCheck(name2))
    
    info = name1 + '-' + name2
    coeff = GetInfo(info, 'bond').split(',')
    
    str1 = "{:>8}{:>6}{:>4}{:>12}{:>13}\n".format(idx_1, idx_2, coeff[0], coeff[1], coeff[2])
    
    return str1
def CountResidue(resName, molList):
    count = 0
    for mol in molList:
        name = mol.getName()
        if name == resName:
            count += 1
    return count

def GetMolNameList(molList):
    nameList = []
    for mol in molList:
        name = mol.getName()
        nameList.append(name)
    return nameList

def ExportBond(pairList, molList, topName):
    output = []
    for i in range(len(pairList)):
        info = GetBondInfo(pairList[i], molList)
        output.append(info)
#        f = open(topName, 'a')
#        f.write(info)
#        f.close()
    return output

def Main(filename, outputName, topName, monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, cutoff, bondCycle):
    df = pd.read_csv(filename, sep='\n', header=None)
    baseList = df.iloc[2:-1].reset_index(drop=True)[0].str.split()
    molName = df.iloc[0][0]
    molNum = df.iloc[1][0]
    molSize = df.iloc[-1][0]
    info = [molName, molNum, molSize]
    
    atomsList = SplitAtom(baseList)
    molList = Atom2Mol(atomsList)
    tmpList = Crosslink(monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, molList, cutoff, bondCycle)
    molList = tmpList[1]
    pairList = tmpList[0]
    print('pairList: ', pairList)
    ExportGRO(info, outputName, molList)
    bondInfo = ExportBond(pairList, molList, topName)
    nameList = GetMolNameList(molList)
#    monNum = CountResidue(monBR, molList)
#    croNum = CountResidue(croBR, molList)
#    moNum = CountResidue(monAR, molList)
#    coNum = CountResidue(croAR, molList)
#    morNum = CountResidue(monAR+'R', molList)
#    corNum = CountResidue(croAR+'R', molList)
    
    
    return nameList, bondInfo
##########################################################
#                         Test                           #
##########################################################
#filename = 'min.gro'
#outputName = 'tst.gro'
#topName = 'tst.top'
#monBR = 'MON'
#monAR = 'MO'
#croBR = 'CRO'
#croAR = 'CO'
#atom1BR = 'C1'
#atom1AR = 'CR'
#atom2BR = 'N1'
#atom2AR = 'NR'
#
#cutoff = 10.
#bondCycle = 2
#
#a=Main(filename, outputName, topName, monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, cutoff, bondCycle)
#print(a)









