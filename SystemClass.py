# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 18:27:02 2018

@author: HuangMing
"""

class SystemInfo(object):
    def __init__(self):
        self.__name__ = ""
        self.__num__ = ""
        self.__size__ = ""
        
    def setName(self, name):
        self.__name__ = name
    
    def setNum(self, num):
        self.__num__ = num
        
    def setSize(self, size):
        self.__size__ = size
    
    def getName(self):
        return self.__name__
    
    def getNum(self):
        return self.__num__
    
    def getSize(self):
        return self.__size__
    
    