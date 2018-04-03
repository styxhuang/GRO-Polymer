# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 14:38:00 2018

@author: huang
"""

class AtomsInfo(object):
    def __init__(self):
        self.__molNum__ = 0
        self.__molName__ = ""
        self.__atomName__ = ""
        self.__globalIndex__ = 0
        self.__localIndex__ = 0
        self.__pos__ = []
        
    def setmolNum(self, num):
        self.__molNum__ = num
        
    def setmolName(self, name):
        self.__molName__ = name
        
    def setatomName(self, name):
        self.__atomName__ = name
    
    def setglobalIndex(self, index):
        self.__globalIndex__ = int(index)
    
    def setlocalIndex(self, index):
        self.__localIndex__ = index
        
    def setPos(self, pos):
        self.__pos__ = pos
        
    
    #################################
    def getmolNum(self):
        return self.__molNum__
    
    def getmolName(self):
        return self.__molName__
    
    def getatomName(self):
        return self.__atomName__
    
    def getglobalIndex(self):
        return self.__globalIndex__
    
    def getlocalIndex(self):
        return self.__localIndex__
    
    def getPos(self):
        return self.__pos__
    
    def outputInfo(self):
        mol = str(self.__molNum__) + self.__molName__
        str1 = '{:>8}{:>7}{:>5}{:>8}{:>8}{:>8}'.format(mol, self.__atomName__, self.__globalIndex__, self.__pos__[0], self.__pos__[1], self.__pos__[2])
        return str1