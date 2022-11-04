import copy
import os
import ROOT as rt

def mkdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory) 

class tnpSample:
    def __init__(self, sampleName, inputPath, outputPath, isMC, lumi=None):
        self.inputPath  = inputPath
        self.outputPath = outputPath
        self.name = sampleName
        self.isMC = isMC
        self.lumi = lumi

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

    def setLumi(self, lumi):
        self.lumi = lumi
        
    def getLumi(self):
        self.lumi

    def clone(self):
        return copy.deepcopy(self)

    def printConfig(self):
        print("="*30)
        print("Name   :", self.getName())
        print("Input  :", self.getInputPath())
        print("Output :", self.getOutputPath())
        print("Lumi   :", self.getLumi())
        print("="*30)

# class tnpVar:
#     def __init__(self, var, hname = None, title = None, xmin = 0, xmax = 0, nbins = -1 ):
#         self.var   = var
#         if title is None :  self.title = var
#         else:               self.title = title
#         self.xmin  = xmin
#         self.xmax  = xmax
#         self.nbins = nbins
#         self.hname = hname
#         self.hist  = None

#     def get_hist(self):
#         if self.nbins > 0:
#             if self.hname is None:  self.hname  = 'h_%' % var
#             self.hist  = rt.TH1F( self.hname, self.title, 
#                                   self.nbins, self.xmin, self.xmax )
#             self.hist.GetXaxis().SetTitle(self.title)
#             self.hist.SetMinimum(0)

#         else:
#             self.hist = None

#         return self.hist

#     def set_hname(self,name):
#         self.hname = name


