import copy
import os
#import ROOT as rt

def mkdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory) 

class tnpSample:
    def __init__(self, sampleName, inputPath, outputPath, isMC):
                 #, lumi=None):
        self.inputPath  = inputPath
        self.outputPath = outputPath
        self.name = sampleName
        self.isMC = isMC
        #self.lumi = lumi

    def isMonteCarlo(self):
        return self.isMC

    def setName(self, name):
        self.name = name

    def getName(self):
        return self.name
        
    def setInputPath(self, inputPath):
        if isinstance(inputPath, str):
            self.inputPath = inputPath
        else:
            print(f"Error in setInputPath() for {self.name}: inputPath type must be string, but it was {type(inputPath)}")

    def getInputPath(self):
        return self.inputPath

    def setOutputPath(self, outputPath):
        if isinstance(outputPath, str):
            self.outputPath = outputPath
        else:
            print(f"Error in setOutputPath() for {self.name}: outputPath type must be string, but it was {type(outputPath)}")

    def getOutputPath(self):
        return self.outputPath

    # def setLumi(self, lumi):
    #     self.lumi = lumi
        
    # def getLumi(self):
    #     self.lumi

    def clone(self):
        return copy.deepcopy(self)

    def printConfig(self):
        print("="*30)
        print("Name   :", self.getName())
        print("Input  :", self.getInputPath())
        print("Output :", self.getOutputPath())
        #print("Lumi   :", self.getLumi())
        print("="*30)
