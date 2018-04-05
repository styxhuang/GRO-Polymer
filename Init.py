# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 09:28:10 2018
Several Hard Code:
    1. Reacted atoms will be changed to from original to C/N
    2. About hydrogen, each bond generation will delete the adjacent H, like C1-H1-H2, after bond creation, it becomes C1-H2
    3. Make difference in one molecule reaction atoms (For reaction, C1 and C2)
    4. If not the 1st atom reacted, then the residue name will add 'R' at the end
    
#TODO: Under the folder Test/gmx-approach/test_0/md has some mismatch warning, needs to be fixed, this is just because the itp file doesn't set properly

@author: HuangMing
"""

import argparse
import fileinput
import subprocess
import os
import sys
from shutil import copyfile
from shutil import move
import pandas as pd
import readGRO

GMX = 'gmx'
def PreMD(monInput, croInput, monNum, croNum, outputName, boxSize): 
    command1 = '{} insert-moleculese -f {} -ci {} -nmol {} -box {} {} {} -o {}'.format(
            GMX, monInput, monInput, int(monNum) - 1, boxSize, boxSize, boxSize, outputName)
    command2 = '{} insert-moleculese -f {} -ci {} -nmol {} -box {} {} {} -o {}'.format(
            GMX, outputName, croInput, croNum, boxSize, boxSize, boxSize, outputName)
    
    subprocess.call(command1, shell=True)
    subprocess.call(command2, shell=True)
    
def MDSimulation():
    command1 = '{} grompp -f em.mdp -c box.gro -p topol.top -o min'
    command2 = '{} mdrun -deffnm min -v'
    command3 = '{} grompp -f nvt.mdp -c min.gro -p topol.top -o nvt'
    command4 = '{} mdrun -deffnm nvt -v'
    subprocess.call(command1, shell=True)
    subprocess.call(command2, shell=True)
    subprocess.call(command3, shell=True)
    subprocess.call(command4, shell=True)

def ExportTop(info, topName, monBR, monAR, croBR, croAR):
#    monBR_num = info[0]
#    monAR_num = info[1]
#    croBR_num = info[2]
#    croAR_num = info[3]
#    bondInfo = info[4]
#    mor_num = info[5]
#    cor_num = info[6]
#    with fileinput.FileInput(topName, inplace=True, backup='.bak') as file:
#        for line in file:
    nameList = info[0]
    bondInfo = info[1]
    
    df= pd.read_csv(topName, sep='\n', header=None)
    idx = df.index[df[0].str.contains('molecules')].tolist() #Assumption, in top file, the molecule column start with monBR
    if len(idx) > 2:
        exit()
    else:
        print('idx: ', idx)
        tmp = df.drop(df.index[int(idx[0])+1:]).reset_index(drop=True)
#        str1 = '{:>5}\t{:>2}'.format(monBR, monBR_num)
#        str2 = '{:>5}\t{:>2}'.format(monAR, monAR_num)
#        str3 = '{:>5}\t{:>2}'.format(croBR, croBR_num)
#        str4 = '{:>5}\t{:>2}'.format(croAR, croAR_num)
#        str5 = '{:>5}\t{:>2}'.format(monAR+'R', mor_num)
#        str6 = '{:>5}\t{:>2}'.format(croAR+'R', cor_num)
        str7 = '\n'
        str8 = '[ intermolecular_interactions ]'
        str9 = '  [ bonds ]'
        
#        Input = [str1, str2, str3, str4, str5, str6, str7, str8, str9]
        for i in range(len(nameList)):
            index1 = idx[0] + i + 1
            str1 = '{}\t{}'.format(nameList[i], 1)
            tmp.loc[index1] = str1
            
        tmp.loc[index1 + 1] = str7
        tmp.loc[index1 + 2] = str8
        tmp.loc[index1 + 3] = str9
        
        for ii in range (len(bondInfo)):
            index2 = index1 + 3 + ii + 1
            print('index', index2)
            str2 = bondInfo[ii]
            tmp.loc[index2] = str2

    f = open('test.top', 'w')        
    for i in range(len(tmp)):
        line = tmp.iloc[i][0] + '\n'
        if '[ atoms ]' in line:
            line = '\n' + line
        f.write(line)
    f.close()
    
#if __name__ == "__main__":
def main():
    parser = argparse.ArgumentParser(description = 'Start Crosslink process')
    
    parser.add_argument('-m', '-mono', action = 'store', dest='mono',
                        help = 'input the monomer mol2 file', metavar = '')
    parser.add_argument('-mL', '-monoLen', action = 'store', dest='monoLen',
                        help = 'input the monomer atoms number w/o H', metavar = '')
    parser.add_argument('-c', '-cross', action = 'store', dest='cross',
                    help = 'input the crosslinker mol2 file', metavar = '')
    parser.add_argument('-cL', '-crossLen', action = 'store', dest='crossLen',
                    help = 'input the crosslinker atoms number w/o H', metavar = '')
    
    parser.add_argument('-b', '-box', dest='box',
                        default = 5, type = int,
                    help = 'input the box size(default box size is 5nm)', metavar = '')
    parser.add_argument('-Nm', '-numberM', dest='num1',
                        default = 4, type = int,
                    help = 'input the monomer number(default monomer number is 4)', metavar = '')
    parser.add_argument('-Nc', '-numberC', dest='num2',
                        default = 2, type = int,
                    help = 'input the crosslinker number(default crosslinker number is 4)', metavar = '')
    parser.add_argument('-r1', '-reactM', action = 'append', dest='r1',
                        default = [], nargs = 2,
                    help = 'input the reactive atoms index in monomer, w/o hydrogen', metavar = '')
    parser.add_argument('-r2', '-reactC', action = 'append', dest='r2',
                        default = [], nargs = 2,
                    help = 'input the reactive atoms index in crosslinker w/o hydrogen', metavar = '')
    parser.add_argument('-cf', '-cutoff', dest='cutoff',
                        default = 5., type = float,
                    help = 'input the radius cutoff(default cutoff equals to 1nm)', metavar = '')
    parser.add_argument('-B', '-Bond', dest='bond',
                        default = 1, type = int,
                    help = 'input the bonds generation each cycle(default bond equals to 1)', metavar = '')
    parser.add_argument('-BT', '-bondTotal', dest='Bonds',
                        default = 1, type = int,
                    help = 'input the total bonds will be generated(default total bonds equals to 1', metavar = '')
       
    args = parser.parse_args()

#    parser.print_help() #For test only
    
    print('Monomer name \t\t\t= ', args.mono)
    print('Monomer atoms number(w/o H) \t= ', args.monoLen)
    print('Crosslinker name \t\t= ', args.cross)
    print('Crosslinker atoms number(w/o H)\t= ', args.crossLen)

    print('Box size \t\t\t= ', args.box)
    print('Monomer number \t\t\t= ', args.num1)
    print('Crosslinker number \t\t= ', args.num2)
    print('Reactive atoms on monomer \t= ', args.r1)
    print('Reactive atoms on crosslinker \t= ', args.r2)
    print('Cutoff \t\t\t\t= ', args.cutoff)
    print('Bonds generation each cycle \t= ', args.bond)
    print('Total bonds will be generated \t= ', args.Bonds)
    
    ###########################################################################
    #               Set Parameters
    ###########################################################################
    monInput = args.mono
    croInput = args.cross    
    monNum = args.num1
    croNum = args.num2

    boxSize = args.box
    initName = 'init.gro'
    PreMD(monInput, croInput, monNum, croNum, initName, boxSize) #setup system

    cutoff = args.cutoff
    bondCycle = args.bond
    totalBonds = args.Bonds
    
    filename = initName
    outputName = 'box.gro'
    topName = 'topolName'
    monBR = 'MON'
    monBR = 'MON'
    monAR = 'MO'
    croBR = 'CRO'
    croAR = 'CO'
    atom1BR = 'C1'
    atom1AR = 'CR'
    atom2BR = 'N1'
    atom2AR = 'NR'
    
    bonds = 0
    i = 0
    
    while bonds < totalBonds:
        i += 1
        command1 = 'cp *itp step{}/'.format(i)
        command2 = 'cp box.gro step{}/'.format(i)
        command3 = 'cp em.mdp step{}/'.format(i)
        command4 = 'cp nvt.mdp step{}/'.format(i)
        command5 = 'cp topol.top step{}/'.format(i)        
        command6 = 'cp nvt.gro ../box.gro'

        os.mkdir('step{}'.format(i))

        subprocess.call(command1, shell=True)
        subprocess.call(command2, shell=True)
        subprocess.call(command3, shell=True)
        subprocess.call(command4, shell=True)
        subprocess.call(command5, shell=True)
        
        os.chdir('step{}'.format(i))
        
        num = readGRO.Main(filename, outputName, topName, monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, cutoff, bondCycle)
        ExportTop(num, topFName, monBR, monAR, croBR, croAR)

        MDSimulation()
        subprocess.call(command6, shell=True)
        bonds += int(bondCycle)

###############################################################################
filename = 'min.gro'
outputName = 'tst.gro'
topName = 'tst.top'
topFName = 'topol.top'
monBR = 'MON'
monAR = 'MO'
croBR = 'CRO'
croAR = 'CO'
atom1BR = 'C1'
atom1AR = 'CR'
atom2BR = 'N1'
atom2AR = 'NR'

cutoff = 10.
bondCycle = 2
info = readGRO.Main(filename, outputName, topName, monBR, monAR, croBR, croAR, atom1BR, atom1AR, atom2BR, atom2AR, cutoff, bondCycle)
ExportTop(info, topFName, monBR, monAR, croBR, croAR)
