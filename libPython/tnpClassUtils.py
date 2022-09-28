import copy
import os
import ROOT as rt

def mkdir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory) 

class tnpSample:
    def __init__(self, sName, path, cut = None, lumi = -1, nEvts = -1, mcTruth = False, puTree = None, isMC = False ):
        self.path = []
        self.name = sName
        self.path.append(path)
        self.cut     = cut
        self.lumi    = lumi
        self.nEvts   = nEvts
        self.mcTruth = mcTruth
        self.puTree  = puTree
        self.isMC    = isMC
        self.weight  = None
        self.tnpTree = None
        self.maxWeight = 1e9

    def set_weight(self,weight):
        self.weight = weight

    def set_maxWeight(self, maxi):
        self.maxWeight = maxi
    
    def set_tnpTree( self, treename):
        self.tnpTree = treename

    def set_puTree(self,puTree):
        self.puTree = puTree

    def set_cut(self,cut):
        self.cut = cut
    
    def set_mcTruth(self,truth = True):
        self.mcTruth = truth

    def dump(self):
        print('**** name: {n} '.format(n=self.name))
        print('  path    : ', self.path)
        print('  tnpTree : ', self.tnpTree)
        if self.isMC:
            print('   --- MC sample --- ')
            print('  nEvts    : ', self.nEvts)
            print('  mcTruth  : ', self.mcTruth)
            print('  puTree   : ', self.puTree)
            print('  weight   : ', self.weight)
        else  :
            print('   --- Data sample --- ')
            print('  lumi     : ', self.lumi)

    def rename(self, newname):
        self.name = newname
    
    def clone(self):
        return copy.deepcopy(self)


    def add_sample(self, sample):
        if self.lumi >= 0  :
            self.lumi = self.lumi + sample.lumi
        if self.nEvts >= 0 :
            self.nEvts = self.nEvts + sample.nEvts
        self.path.extend( sample.path )


class tnpSampleNew:
    def __init__(self, sampleName, inputPath, outputPath, isMC, lumi):
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
        if self.isMC:
            print("Lumi   :", self.getLumi())
        print("="*30)

class tnpVar:
    def __init__(self, var, hname = None, title = None, xmin = 0, xmax = 0, nbins = -1 ):
        self.var   = var
        if title is None :  self.title = var
        else:               self.title = title
        self.xmin  = xmin
        self.xmax  = xmax
        self.nbins = nbins
        self.hname = hname
        self.hist  = None

    def get_hist(self):
        if self.nbins > 0:
            if self.hname is None:  self.hname  = 'h_%' % var
            self.hist  = rt.TH1F( self.hname, self.title, 
                                  self.nbins, self.xmin, self.xmax )
            self.hist.GetXaxis().SetTitle(self.title)
            self.hist.SetMinimum(0)

        else:
            self.hist = None

        return self.hist

    def set_hname(self,name):
        self.hname = name


