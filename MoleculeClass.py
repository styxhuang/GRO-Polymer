# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 12:53:50 2018

@author: huang
"""

class MoleculeInfo(object):
    def __init__(self):
        self.__index__ = 0
        self.__name__ = ""
        self.__atoms__ = []
        
    def setIndex(self, index):
        self.__index__ = index
    
    def setName(self, name):
        self.__name__ = name
        
    def setAtoms(self, atom):
        self.__atoms__.append(atom)
        
    def getIndex(self):
        return self.__index__
    
    def getName(self):
        return self.__name__
    
    def getAtoms(self):
        return self.__atoms__
    
    def outputInfo(self):
        print("Index: ", self.__index__, "\tName: ", self.__name__, "\t", len(self.__atoms__), "\t", "Reactive atoms: ", self.__reactAtoms__[0])
        for i in range (len(self.__atoms__)):
            index = self.__atoms__[i].getIndex()
            name = self.__atoms__[i].getAtomName()
            subName = self.__atoms__[i].getSubName()
            print(index, "\t", name, "\t", subName)
            
            