import ROOT

#for so in ['histFitter_C.so', 'histScaleFitter_C.so', 'RooCBExGaussShape_cc.so', 'RooCMSShape_cc.so']:
#
#    if so not in ROOT.gSystem.GetLibraries():
#        print('need to build', so)
#        ROOT.gROOT.ProcessLine('.L ./libCpp/{soC}+'.format(soC=so.replace('_C.so','.C').replace('_cc.so','.cc')))

##ROOT.gROOT.ProcessLine('.L ./libCpp/histFitter.C+')
#ROOT.gROOT.LoadMacro('./libCpp/histFitter.C+')
##ROOT.gROOT.LoadMacro('./libCpp/histScaleFitter.C+')
#ROOT.gROOT.LoadMacro('./libCpp/RooCBExGaussShape.cc+')
#ROOT.gROOT.LoadMacro('./libCpp/RooCMSShape.cc+')
#ROOT.gROOT.SetBatch(1)

from ROOT import RooFit,RooFitResult
from ROOT import tnpFitter,scaleFitter

import re
import math


minPtForSwitch = 70

def ptMin( tnpBin ):
    ptmin = 1
    if tnpBin['name'].find('pt_') >= 0:
        ptmin = float(tnpBin['name'].split('pt_')[1].split('p')[0])
    elif tnpBin['name'].find('et_') >= 0:
        ptmin = float(tnpBin['name'].split('et_')[1].split('p')[0])
    return ptmin

def createWorkspaceForAltSig( sample, tnpBin, tnpWorkspaceParam, refResample=-1 ):

    if sample.isMonteCarlo():
        return tnpWorkspaceParam
    
    fileref = sample.mcRef.altSigFit
    filemc  = ROOT.TFile(fileref,'read')

    if refResample<0:
        fitresP = filemc.Get( '%s_resP' % tnpBin['name']  )
    else:
        fitresP = filemc.Get( '%s_resP_Stat%d' % (tnpBin['name'],refResample)  )
    fitresF = filemc.Get( '%s_resF' % tnpBin['name'] )

    listOfParamP = ['nP', 'alphaP', 'sigmaP', 'sigmaP_2']
    listOfParamF = ['nF', 'alphaF', 'sigmaF', 'sigmaF_2']

    # set central value of signal parameters as in MC alt sig fit, but do not fix them
    # while those for background from the nominal fit in data (but only for failing probes, passing ones are good)
    
    # first for failing probes and signal parameters
    fitPar = fitresF.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        #print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        x = re.compile('%s.*?' % pName)
        for par in listOfParamF:
            if pName == par:
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                ir = listToRM[0] # should always be only 1 element, otherwise adding it back below becomes a problem
                #print(f">>>>> old {ir}")
                tnpWorkspaceParam.remove(ir)
                parRange = ir.split("[")[1].split("]")[0].split(",")  # get elements within square brackets, "ir" is like "name[value,low,high]"
                parRange = ",".join(parRange[1:]) # concatenate everything except the first element
                new_ir = "%s[%2.3f,%s]" % (pName, fitPar[ipar].getVal(), parRange)
                #print(f">>>>> new {new_ir}")
                tnpWorkspaceParam.append( new_ir )
                    
    # now for passing probes
    fitPar = fitresP.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        #print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        x = re.compile('%s.*?' % pName)
        for par in listOfParamP:
            if pName == par:
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                ir = listToRM[0] # should always be only 1 element, otherwise adding it back below becomes a problem
                #print(f">>>>> old {ir}")
                tnpWorkspaceParam.remove(ir)
                parRange = ir.split("[")[1].split("]")[0].split(",")  # get elements within square brackets, "ir" is like "name[value,low,high]"
                parRange = ",".join(parRange[1:]) # concatenate everything except the first element
                new_ir = "%s[%2.3f,%s]" % (pName, fitPar[ipar].getVal(), parRange)
                #print(f">>>>> new {new_ir}")
                tnpWorkspaceParam.append( new_ir )
                
    filemc.Close()

    # print(">>>>>")
    # print(sample)
    # print(sample.nominalFit)
    # print(">>>>>")
    filerefData = sample.nominalFit
    filedata  = ROOT.TFile(filerefData,'read')
    fitresF = filedata.Get( '%s_resF' % tnpBin['name'] )
    listOfBkgParamF = ['acmsF', 'betaF', 'gammaF']
    
    # failing probes and background parameters
    fitPar = fitresF.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        #print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        x = re.compile('%s.*?' % pName)
        for par in listOfBkgParamF:
            if pName == par:
                listToRM = list(filter(x.match, tnpWorkspaceParam))
                ir = listToRM[0] # should always be only 1 element, otherwise adding it back below becomes a problem
                #print(f">>>>> old {ir}")
                tnpWorkspaceParam.remove(ir)
                parRange = ir.split("[")[1].split("]")[0].split(",")  # get elements within square brackets, "ir" is like "name[value,low,high]"
                parRange = ",".join(parRange[1:]) # concatenate everything except the first element
                new_ir = "%s[%2.3f,%s]" % (pName, fitPar[ipar].getVal(), parRange)
                #print(f">>>>> new {new_ir}")
                tnpWorkspaceParam.append( new_ir )
                                    
    filedata.Close()

    
    return tnpWorkspaceParam


#############################################################
########## nominal fitter
#############################################################
def histFitterNominal( sample, tnpBin, tnpWorkspaceParam, massbins=60, massmin=60, massmax=120, useAllTemplateForFail=False, maxFailIntegralToUseAllProbe=-1):
        
    print('------- now nominal fitting bin:')
    for i in tnpBin:
        print(i, tnpBin[i])
    tnpWorkspaceFunc = [
        "Gaussian::sigResPass(x,meanP,sigmaP)",
        "Gaussian::sigResFail(x,meanF,sigmaF)",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)
    
    ## init fitter
    infile = ROOT.TFile( sample.getOutputPath(), "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'], massbins, massmin, massmax )
    infile.Close()

    ## setup
    ## make configurable from outside
    fitter.useMinos()
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    fitter.setOutputFile( sample.nominalFit+'_bin_'+tnpBin['name'])
    
    ## generated Z LineShape
    ## for high pT change the failing spectra to any probe to get statistics
    fileTruth  = ROOT.TFile(sample.mcRef.getOutputPath(),'read')
    histZLineShapeP = fileTruth.Get('%s_Pass'%tnpBin['name'])
    histZLineShapeF = fileTruth.Get('%s_Fail'%tnpBin['name'])
        
    if useAllTemplateForFail or ptMin(tnpBin) > minPtForSwitch:
        if maxFailIntegralToUseAllProbe > 0 and histZLineShapeF.Integral() < maxFailIntegralToUseAllProbe:
            histZLineShapeF.Add(histZLineShapeP)
    
    fitter.setZLineShapes(histZLineShapeP,histZLineShapeF)
    fileTruth.Close()

    # python3rootfile.cd()
    ### set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace, sample.isMonteCarlo() )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_eta','#eta')
    title = title.replace('probe_pt','p_{T}')
    fitter.fits(title)
    # python3rootfile.Close()

#############################################################
########## alternate signal fitter
#############################################################
def histFitterAltSig( sample, tnpBin, tnpWorkspaceParam, massbins=60, massmin=60, massmax=120, altSignalFail=False):

    tnpWorkspacePar = createWorkspaceForAltSig( sample,  tnpBin, tnpWorkspaceParam )

    ### tricky: use n < 0 for high pT bin (so need to remove param and add it back)
    ptmin = ptMin(tnpBin)        
    
    # FIXME: is the tail not overridden inside createWorkspaceForAltSig?
    if altSignalFail:
        tnpWorkspaceFunc = [
            "tailLeft[%d]" % (-1 if ptmin >= 35 else 1),
            "RooCBExGaussShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)",
            "Gaussian::sigResFail(x,meanF,sigmaF)",
            "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]
    else:
        tnpWorkspaceFunc = [
            "tailLeft[%d]" % (-1 if ptmin >= 35 else 1),
            "RooCBExGaussShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)",
            "RooCBExGaussShape::sigResFail(x,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),alphaF,nF, expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),tailLeft)",
            "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
            "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]


    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspacePar)
    tnpWorkspace.extend(tnpWorkspaceFunc)
        
    ## init fitter
    infile = ROOT.TFile( sample.getOutputPath(), "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    ## for high pT change the failing spectra to passing probe to get statistics 
    ## MC only: this is to get MC parameters in data fit!
    if sample.isMonteCarlo() and ptMin( tnpBin ) > minPtForSwitch:     
        hF = infile.Get('%s_Pass' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'], massbins, massmin, massmax )
#    fitter.fixSigmaFtoSigmaP()
    infile.Close()

    ## setup
    ## make configurable from outside
    fitter.useMinos()
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    fitter.setOutputFile( sample.altSigFit+'_bin_'+tnpBin['name'])
    
    ## generated Z LineShape
    fileTruth = ROOT.TFile('etc/inputs/ZeeGenLevel.root','read')
    histZLineShape = fileTruth.Get('Mass')
    fitter.setZLineShapes(histZLineShape,histZLineShape)
    fileTruth.Close()

    ### set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace, sample.isMonteCarlo() )

    title = tnpBin['title'].replace(';',' - ')
    fitter.fits(title)

    #rootfile.Close()


#############################################################
########## alternate background fitter
#############################################################
def histFitterAltBkg( sample, tnpBin, tnpWorkspaceParam, massbins=60, massmin=60, massmax=120, useAllTemplateForFail=False, maxFailIntegralToUseAllProbe=-1):

    tnpWorkspaceFunc = [
        "Gaussian::sigResPass(x,meanP,sigmaP)",
        "Gaussian::sigResFail(x,meanF,sigmaF)",
        "Exponential::bkgPass(x, alphaP)",
        "Exponential::bkgFail(x, alphaF)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)
            
    ## init fitter
    infile = ROOT.TFile(sample.getOutputPath(),'read')
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'], massbins, massmin, massmax )
    infile.Close()

    ## setup
    ## make configurable from outside
    fitter.useMinos()
    fitter.setPassStrategy(2)
    fitter.setFailStrategy(2)
    fitter.setPrintLevel(-1)
    fitter.setOutputFile(sample.altBkgFit+'_bin_'+tnpBin['name'])

    ## generated Z LineShape
    ## for high pT change the failing spectra to any probe to get statistics
    fileTruth = ROOT.TFile(sample.mcRef.getOutputPath(),'read')
    histZLineShapeP = fileTruth.Get('%s_Pass'%tnpBin['name'])
    histZLineShapeF = fileTruth.Get('%s_Fail'%tnpBin['name'])
    if useAllTemplateForFail or ptMin(tnpBin) > minPtForSwitch:
        if maxFailIntegralToUseAllProbe > 0 and histZLineShapeF.Integral() < maxFailIntegralToUseAllProbe:
            histZLineShapeF.Add(histZLineShapeP)

    fitter.setZLineShapes(histZLineShapeP,histZLineShapeF)
    fileTruth.Close()

    ### set workspace
    workspace = ROOT.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace, sample.isMonteCarlo() )

    title = tnpBin['title'].replace(';',' - ')
    fitter.fits(title)
    #rootfile.Close()


