import ROOT as rt

#for so in ['histFitter_C.so', 'histScaleFitter_C.so', 'RooCBExGaussShape_cc.so', 'RooCMSShape_cc.so']:
#
#    if so not in rt.gSystem.GetLibraries():
#        print('need to build', so)
#        rt.gROOT.ProcessLine('.L ./libCpp/{soC}+'.format(soC=so.replace('_C.so','.C').replace('_cc.so','.cc')))

#rt.gROOT.ProcessLine('.L ./libCpp/histFitter.C+')
rt.gROOT.LoadMacro('./libCpp/histFitter.C+')
rt.gROOT.LoadMacro('./libCpp/histScaleFitter.C+')
rt.gROOT.LoadMacro('./libCpp/RooCBExGaussShape.cc+')
rt.gROOT.LoadMacro('./libCpp/RooCMSShape.cc+')
rt.gROOT.SetBatch(1)

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

def createWorkspaceForAltSig( sample, tnpBin, tnpWorkspaceParam, tnpFit=True, refResample=-1 ):

    ### tricky: use n < 0 for high pT bin (so need to remove param and add it back)
    cbNList = ['tailLeft']
    ptmin = ptMin(tnpBin)        
    if ptmin >= 35 :
        for par in cbNList:
            for ip in range(len(tnpWorkspaceParam)):
                x=re.compile('%s.*?' % par)
                listToRM = filter(x.match, tnpWorkspaceParam)
                for ir in listToRM :
                    print('**** remove', ir)
                    tnpWorkspaceParam.remove(ir)                    
            tnpWorkspaceParam.append( 'tailLeft[-1]' )

    if sample.isMC:
        return tnpWorkspaceParam

    
    fileref = sample.mcRef.altSigFit
    filemc  = rt.TFile(fileref,'read')

    from ROOT import RooFit,RooFitResult
    if refResample<0:
        fitresP = filemc.Get( '%s_resP' % tnpBin['name']  )
    else:
        fitresP = filemc.Get( '%s_resP_Stat%d' % (tnpBin['name'],refResample)  )
    if tnpFit: fitresF = filemc.Get( '%s_resF' % tnpBin['name'] )

    listOfParam = ['nF','alphaF','nP','alphaP','sigmaP','sigmaF','sigmaP_2','sigmaF_2']
    
    if tnpFit:
        fitPar = fitresF.floatParsFinal()
        for ipar in range(len(fitPar)):
            pName = fitPar[ipar].GetName()
            print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
            for par in listOfParam:
                if pName == par:
                    x=re.compile('%s.*?' % pName)
                    listToRM = filter(x.match, tnpWorkspaceParam)
                    for ir in listToRM :
                        tnpWorkspaceParam.remove(ir)                    
                    tnpWorkspaceParam.append( '%s[%2.3f]' % (pName,fitPar[ipar].getVal()) )
                              
  
    fitPar = fitresP.floatParsFinal()
    for ipar in range(len(fitPar)):
        pName = fitPar[ipar].GetName()
        print('{n}[{f:.3f}]'.format(n=pName,f=fitPar[ipar].getVal()))
        for par in listOfParam:
            if pName == par:
                x=re.compile('%s.*?' % pName)
                listToRM = filter(x.match, tnpWorkspaceParam)
                for ir in listToRM :
                    tnpWorkspaceParam.remove(ir)
                tnpWorkspaceParam.append( '%s[%2.3f]' % (pName,fitPar[ipar].getVal()) )

    filemc.Close()

    return tnpWorkspaceParam


#############################################################
########## nominal fitter
#############################################################
def histFitterNominal( sample, tnpBin, tnpWorkspaceParam ):
        
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
    infile = rt.TFile( sample.histFile, "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'] )
    infile.Close()

    ## setup
    fitter.useMinos()
    # python3rootfile = rt.TFile(sample.nominalFit+'_bin_'+tnpBin['name'],'update')
    fitter.setOutputFile( sample.nominalFit+'_bin_'+tnpBin['name'])
    
    ## generated Z LineShape
    ## for high pT change the failing spectra to any probe to get statistics
    fileTruth  = rt.TFile(sample.mcRef.histFile,'read')
    histZLineShapeP = fileTruth.Get('%s_Pass'%tnpBin['name'])
    histZLineShapeF = fileTruth.Get('%s_Fail'%tnpBin['name'])
    if ptMin( tnpBin ) > minPtForSwitch: 
        histZLineShapeF = fileTruth.Get('%s_Pass'%tnpBin['name'])
#        fitter.fixSigmaFtoSigmaP()
    fitter.setZLineShapes(histZLineShapeP,histZLineShapeF)

    fileTruth.Close()

    # python3rootfile.cd()
    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)
    # python3rootfile.Close()

#############################################################
########## nominal scale fitter
#############################################################
def histScaleFitterNominal( sample, tnpBin, tnpWorkspaceParam, resample, batch=False ):
        
    tnpWorkspaceFunc = [
        "Gaussian::sigResPass(x,meanP,sigmaP)",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspaceParam)
    tnpWorkspace.extend(tnpWorkspaceFunc)
    
    ## init fitter
    infile = rt.TFile( getattr(sample, 'histFile%d' % resample), "read")
    hP = infile.Get('%s_Stat%d' % (tnpBin['name'], resample) )
    fitter = scaleFitter( hP, tnpBin['name'], resample )
    infile.Close()

    ## setup
    fitter.useMinos()
    rootfilename = sample.nominalFit if not batch else sample.nominalFit.replace('.root','_Stat%d.root' % resample)
    rootfile = rt.TFile(rootfilename,'update')
    fitter.setOutputFile( rootfile )
    
    ## generated Z LineShape
    fileTruth  = rt.TFile(getattr(sample.mcRef,'histFile%d' % resample),'read')
    histZLineShapeP = fileTruth.Get('%s_Stat%d' % (tnpBin['name'],resample))
    fitter.setZLineShape(histZLineShapeP)

    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)
    rootfile.Close()


#############################################################
########## alternate signal fitter
#############################################################
def histFitterAltSig( sample, tnpBin, tnpWorkspaceParam ):

    tnpWorkspacePar = createWorkspaceForAltSig( sample,  tnpBin, tnpWorkspaceParam )

    tnpWorkspaceFunc = [
        "tailLeft[1]",
        "RooCBExGaussShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)",
        "RooCBExGaussShape::sigResFail(x,meanF,expr('sqrt(sigmaF*sigmaF+sosF*sosF)',{sigmaF,sosF}),alphaF,nF, expr('sqrt(sigmaF_2*sigmaF_2+sosF*sosF)',{sigmaF_2,sosF}),tailLeft)",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        "RooCMSShape::bkgFail(x, acmsF, betaF, gammaF, peakF)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspacePar)
    tnpWorkspace.extend(tnpWorkspaceFunc)
        
    ## init fitter
    infile = rt.TFile( sample.histFile, "read")
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    ## for high pT change the failing spectra to passing probe to get statistics 
    ## MC only: this is to get MC parameters in data fit!
    if sample.isMC and ptMin( tnpBin ) > minPtForSwitch:     
        hF = infile.Get('%s_Pass' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'] )
#    fitter.fixSigmaFtoSigmaP()
    infile.Close()

    ## setup
    ## parallel rootfile = rt.TFile(sample.altSigFit,'update')
    #python3rootfile = rt.TFile(sample.altSigFit+'_bin_'+tnpBin['name'],'update')
    #python3fitter.setOutputFile( rootfile )
    fitter.setOutputFile( sample.altSigFit+'_bin_'+tnpBin['name'])
    
    ## generated Z LineShape
    fileTruth = rt.TFile('etc/inputs/ZeeGenLevel.root','read')
    histZLineShape = fileTruth.Get('Mass')
    fitter.setZLineShapes(histZLineShape,histZLineShape)
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)

    #rootfile.Close()


#############################################################
########## alternate signal scale fitter
#############################################################
def histScaleFitterAltSig( sample, tnpBin, tnpWorkspaceParam, resample, batch=False ):

    # here a random MC replica is used just as nominal
    tnpWorkspacePar = createWorkspaceForAltSig( sample,  tnpBin, tnpWorkspaceParam, tnpFit=False, refResample=resample )

    tnpWorkspaceFunc = [
        "tailLeft[1]",
        "RooCBExGaussShape::sigResPass(x,meanP,expr('sqrt(sigmaP*sigmaP+sosP*sosP)',{sigmaP,sosP}),alphaP,nP, expr('sqrt(sigmaP_2*sigmaP_2+sosP*sosP)',{sigmaP_2,sosP}),tailLeft)",
        "RooCMSShape::bkgPass(x, acmsP, betaP, gammaP, peakP)",
        ]

    tnpWorkspace = []
    tnpWorkspace.extend(tnpWorkspacePar)
    tnpWorkspace.extend(tnpWorkspaceFunc)
        
    ## init fitter
    infile = rt.TFile( getattr(sample, 'histFile%d' % resample), "read")
    print("Getting the histograms to fit the Z lineshape from: ",infile)
    hP = infile.Get('%s_Stat%d' % (tnpBin['name'],resample))
    fitter = scaleFitter( hP, tnpBin['name'], resample )
    infile.Close()

    ## setup
    if batch: sample.altSigFit = sample.altSigFit.replace('.root','_%s.root' % tnpBin['name'])
    #python3rootfile = rt.TFile(sample.altSigFit,'update')
    fitter.setOutputFile( sample.altSigFit+'_bin_'+tnpBin['name'])
    
    ## generated Z LineShape
    fileTruth = rt.TFile('etc/inputs/ZeeGenLevel.root','read')
    histZLineShape = fileTruth.Get('Mass')
    fitter.setZLineShape(histZLineShape)
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)

    #python3rootfile.Close()



#############################################################
########## alternate background fitter
#############################################################
def histFitterAltBkg( sample, tnpBin, tnpWorkspaceParam ):

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
    infile = rt.TFile(sample.histFile,'read')
    hP = infile.Get('%s_Pass' % tnpBin['name'] )
    hF = infile.Get('%s_Fail' % tnpBin['name'] )
    fitter = tnpFitter( hP, hF, tnpBin['name'] )
    infile.Close()

    ## setup
    ## parallel rootfile = rt.TFile(sample.altBkgFit,'update')
    rootfile = rt.TFile(sample.altBkgFit+'_bin_'+tnpBin['name'],'update')
    fitter.setOutputFile( rootfile )
#    fitter.setFitRange(65,115)

    ## generated Z LineShape
    ## for high pT change the failing spectra to any probe to get statistics
    fileTruth = rt.TFile(sample.mcRef.histFile,'read')
    histZLineShapeP = fileTruth.Get('%s_Pass'%tnpBin['name'])
    histZLineShapeF = fileTruth.Get('%s_Fail'%tnpBin['name'])
    if ptMin( tnpBin ) > minPtForSwitch: 
        histZLineShapeF = fileTruth.Get('%s_Pass'%tnpBin['name'])
#        fitter.fixSigmaFtoSigmaP()
    fitter.setZLineShapes(histZLineShapeP,histZLineShapeF)
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)
    rootfile.Close()


#############################################################
########## alternate background scale fitter
#############################################################
def histScaleFitterAltBkg( sample, tnpBin, tnpWorkspaceParam, resample ):

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
    infile = rt.TFile(getattr(sample, 'histFile%d' % resample),'read')
    hP = infile.Get('%s_Stat%d' % (tnpBin['name'],resample))
    fitter = scaleFitter( hP, tnpBin['name'], resample )
    infile.Close()

    ## setup
    rootfile = rt.TFile(sample.altBkgFit,'update')
    fitter.setOutputFile( rootfile )
#    fitter.setFitRange(65,115)

    ## generated Z LineShape
    fileTruth = rt.TFile(getattr(sample.mcRef,'histFile%d' % resample),'read')
    histZLineShapeP = fileTruth.Get('%s_Stat%d' % (tnpBin['name'],resample))
    fitter.setZLineShape(histZLineShapeP)
    fileTruth.Close()

    ### set workspace
    workspace = rt.vector("string")()
    for iw in tnpWorkspace:
        workspace.push_back(iw)
    fitter.setWorkspace( workspace )

    title = tnpBin['title'].replace(';',' - ')
    title = title.replace('probe_sc_eta','#eta_{SC}')
    title = title.replace('probe_Ele_pt','p_{T}')
    fitter.fits(sample.mcTruth,title)
    rootfile.Close()


