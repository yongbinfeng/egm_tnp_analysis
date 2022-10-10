import ROOT, math, array, ctypes, copy
from ROOT import RooFit,RooFitResult
import numpy as np
import os.path
from . import fitUtils
# import *
import functools
#from fitSimultaneousUtils import *

ROOT.ROOT.EnableImplicitMT()
#ROOT.ROOT.TTreeProcessorMT.SetMaxTasksPerFilePerWorker(1)

#ROOT.gROOT.LoadMacro('pileupWeights.C+')

def removeNegativeBins(h):
    for i in range(1,h.GetNbinsX()+1):
        if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)


def makePassFailHistograms(sample, bins, bindef, var ):

    probe_binning_eta, probe_binning_pt = bindef['eta']['bins'], bindef['pt']['bins']
    probe_var_eta, probe_var_pt         = bindef['eta']['var'] , bindef['pt']['var']

    probe_binning_pt  = array.array('d', probe_binning_pt)
    probe_binning_eta = array.array('d', probe_binning_eta)

    probe_var_pass_pt, probe_var_pass_eta = probe_var_pt, probe_var_eta

    #binning_mass = array.array('d', [var['min'] + i*(var['max']-var['min'])/var['nbins'] for i in range(var['nbins']+1)])

    #print("sample.getInputPath() = ",sample.getInputPath()) 
    p = sample.getInputPath() 
    #print("p = ",p)
    infile = ROOT.TFile(p)
    #print(infile.ls())
    h_tmp_pass = infile.Get("pass_" + sample.getName())
    h_tmp_fail = infile.Get("fail_" + sample.getName())
    altPass = "pass_" + sample.getName() + "_alt"
    keyNames = [k.GetName() for k in infile.GetListOfKeys()]
    h_tmp_pass_alt = None # will be needed for tracking when using all probes to form failing probe MC template to fit data, in this case the alternate version of passing probe made with  standalone variables is needed 
    if altPass in keyNames:
        h_tmp_pass_alt = infile.Get(altPass)        
        
    outfile = ROOT.TFile(sample.getOutputPath(), 'recreate')
    
    #print('this is the integral of the passing', h_tmp_pass.Integral())
    #print('this is the integral of the failing', h_tmp_fail.Integral())
    
    for ii,ib in enumerate(bins):
        h_name = ib['name' ]
        h_title= ib['title']

        tmp_valpt  = ib['vars'][probe_var_pt ]['min']
        tmp_valeta = ib['vars'][probe_var_eta]['min']

        epsilon = 0.001 # safety thing when picking the bin edges using FindFixBin
        ibin_pt  = h_tmp_pass.GetYaxis().FindFixBin(tmp_valpt + epsilon)
        ibin_eta = h_tmp_pass.GetZaxis().FindFixBin(tmp_valeta+ epsilon)
        #print('i am at ibin_pt  {ipt}  and pt  {pt:.1f} '.format(ipt = ibin_pt , pt = tmp_valpt ))
        #print('i am at ibin_eta {ieta} and eta {eta:.1f}'.format(ieta= ibin_eta, eta= tmp_valeta))

        h_pass = h_tmp_pass.ProjectionX(h_name+'_Pass', ibin_pt, ibin_pt, ibin_eta, ibin_eta)
        h_fail = h_tmp_fail.ProjectionX(h_name+'_Fail', ibin_pt, ibin_pt, ibin_eta, ibin_eta)
        h_pass.SetTitle(h_title+' passing')
        h_fail.SetTitle(h_title+' failing')
        removeNegativeBins(h_pass)
        removeNegativeBins(h_fail)
        h_pass.Write(h_pass.GetName())
        h_fail.Write(h_fail.GetName())
        if h_tmp_pass_alt:
            h_pass_alt = h_tmp_pass_alt.ProjectionX(h_name+'_Pass_alt', ibin_pt, ibin_pt, ibin_eta, ibin_eta)
            h_pass_alt.SetTitle(h_title+' passing alternate')
            removeNegativeBins(h_pass_alt)
            h_pass_alt.Write(h_pass_alt.GetName())
            
    outfile.Close()


def histPlotter( rootfile, tnpBin, plotDir, replica=-1, verbosePlotting=True ):

    if verbosePlotting:
        if replica<0:
            #print('  get canvas: ' , '{c}_Canv'.format(c=tnpBin['name']))
            c = rootfile.Get( '%s_Canv' % tnpBin['name'] )
            c.Print( '%s/%s.png' % (plotDir,tnpBin['name']))
            c.Print( '%s/%s.pdf' % (plotDir,tnpBin['name']))
        else:
            #print('  get canvas: ' , '{c}_Canv_Stat{r}'.format(c=tnpBin['name'],r=replica))
            c = rootfile.Get( '%s_Canv_Stat%d' % (tnpBin['name'],replica) )
            c.Print( '%s/%s_Stat%d.png' % (plotDir,tnpBin['name'],replica))
            c.Print( '%s/%s_Stat%d.pdf' % (plotDir,tnpBin['name'],replica))
        


def computeEffi( n1,n2,e1,e2):
    effout = []
    if (n1+n2):
        eff   = n1/(n1+n2)
        nTot = n1+n2
        e_eff = math.sqrt(e1*e1*n2*n2+e2*e2*n1*n1) / (nTot*nTot)
        #if e_eff < 0.001 : e_eff = 0.001
    else:
        eff, e_eff = 1.1, 0.01

    effout.append(eff)
    effout.append(e_eff)
    
    return effout


def getAllEffi( info, bindef ):
    effis = {}
    if not info['mcNominal'] is None and os.path.isfile(info['mcNominal']):
        rootfile = ROOT.TFile( info['mcNominal'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])
        bin1 = 1
        bin2 = hP.GetXaxis().GetNbins()
        eP = ctypes.c_double(-1.0)
        eF = ctypes.c_double(-1.0)
        nP = hP.IntegralAndError(bin1,bin2,eP)
        nF = hF.IntegralAndError(bin1,bin2,eF)

        effis['mcNominal'] = computeEffi(nP,nF,float(eP.value),float(eF.value)) +[nP, nF, float(eP.value), float(eF.value)]
        rootfile.Close()
    else: effis['mcNominal'] = [-1,-1]

    if not info['mcNominal_fit'] is None and os.path.isfile(info['mcNominal_fit']):
        rootfile = ROOT.TFile( info['mcNominal_fit'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_mcNominal_fit'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name'] ))

        fitP = fitresP.floatParsFinal().find('nSigP')
        fitF = fitresF.floatParsFinal().find('nSigF')
        
        nP = fitP.getVal()
        nF = fitF.getVal()
        eP = fitP.getError()
        eF = fitF.getError()

        effis['mcNominal_fit'] = computeEffi(nP,nF,eP,eF) +[nP, nF, float(eP), float(eF)]
        rootfile.Close()
    else: 
        effis['mcNominal_fit'] = [-1,-1]
        effis['canv_mcNominal_fit'] = None
    

    if not info['tagSel'] is None and os.path.isfile(info['tagSel']):
        rootfile = ROOT.TFile( info['tagSel'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])
        bin1 = 1
        bin2 = hP.GetXaxis().GetNbins()
        eP = ROOT.Double(-1.0)
        eF = ROOT.Double(-1.0)
        nP = hP.IntegralAndError(bin1,bin2,eP)
        nF = hF.IntegralAndError(bin1,bin2,eF)

        effis['tagSel'] = computeEffi(nP,nF,eP,eF)
        rootfile.Close()
    else: effis['tagSel'] = [-1,-1]
        
    if not info['mcAlt'] is None and os.path.isfile(info['mcAlt']):
        rootfile = ROOT.TFile( info['mcAlt'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_mcAlt'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name'] ))

        fitP = fitresP.floatParsFinal().find('nSigP')
        fitF = fitresF.floatParsFinal().find('nSigF')
        
        nP = fitP.getVal()
        nF = fitF.getVal()
        eP = fitP.getError()
        eF = fitF.getError()

        effis['mcAlt'] = computeEffi(nP,nF,eP,eF) +[nP, nF, float(eP), float(eF)]
        rootfile.Close()
    else: 
        effis['mcAlt'] = [-1,-1]
        effis['canv_mcAlt'] = None

    
    if not info['mcAltBkg'] is None and os.path.isfile(info['mcAltBkg']):
        rootfile = ROOT.TFile( info['mcAltBkg'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_mcAltBkg'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name'] ))

        fitP = fitresP.floatParsFinal().find('nSigP')
        fitF = fitresF.floatParsFinal().find('nSigF')
        
        nP = fitP.getVal()
        nF = fitF.getVal()
        eP = fitP.getError()
        eF = fitF.getError()

        effis['mcAltBkg'] = computeEffi(nP,nF,eP,eF) +[nP, nF, float(eP), float(eF)]
        rootfile.Close()
    else: 
        effis['mcAltBkg'] = [-1,-1]
        effis['canv_mcAltBkg'] = None

    if not info['dataNominal'] is None and os.path.isfile(info['dataNominal']) :
        rootfile = ROOT.TFile( info['dataNominal'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_dataNominal'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name']  ))

        fitP = fitresP.floatParsFinal().find('nSigP')
        fitF = fitresF.floatParsFinal().find('nSigF')
        
        nP = fitP.getVal()
        nF = fitF.getVal()
        eP = fitP.getError()
        eF = fitF.getError()
        rootfile.Close()

        rootfile = ROOT.TFile( info['data'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])

        if eP > math.sqrt(hP.Integral()) : eP = math.sqrt(hP.Integral())
        if eF > math.sqrt(hF.Integral()) : eF = math.sqrt(hF.Integral())
        rootfile.Close()

        effis['dataNominal'] = computeEffi(nP,nF,eP,eF) +[nP, nF, eP, eF]
    else:
        effis['dataNominal'] = [-1,-1]
        effis['canv_dataNominal'] = None

    if not info['dataAltSig'] is None and os.path.isfile(info['dataAltSig']) :
        rootfile = ROOT.TFile( info['dataAltSig'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_dataAltSig'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name'] ))

        nP = fitresP.floatParsFinal().find('nSigP').getVal()
        nF = fitresF.floatParsFinal().find('nSigF').getVal()
        eP = fitresP.floatParsFinal().find('nSigP').getError()
        eF = fitresF.floatParsFinal().find('nSigF').getError()
        rootfile.Close()

        rootfile = ROOT.TFile( info['data'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])

        if eP > math.sqrt(hP.Integral()) : eP = math.sqrt(hP.Integral())
        if eF > math.sqrt(hF.Integral()) : eF = math.sqrt(hF.Integral())
        rootfile.Close()

        effis['dataAltSig'] = computeEffi(nP,nF,eP,eF) +[nP, nF, eP, eF]

    else:
        effis['dataAltSig'] = [-1,-1]
        effis['canv_dataAltSig'] = None

    if not info['dataAltBkg'] is None and os.path.isfile(info['dataAltBkg']):
        rootfile = ROOT.TFile( info['dataAltBkg'], 'read' )
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )
        effis['canv_dataAltBkg'] = copy.deepcopy(rootfile.Get( '%s_Canv' % bindef['name'] ))

        nP = fitresP.floatParsFinal().find('nSigP').getVal()
        nF = fitresF.floatParsFinal().find('nSigF').getVal()
        eP = fitresP.floatParsFinal().find('nSigP').getError()
        eF = fitresF.floatParsFinal().find('nSigF').getError()
        rootfile.Close()

        rootfile = ROOT.TFile( info['data'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])

        if eP > math.sqrt(hP.Integral()) : eP = math.sqrt(hP.Integral())
        if eF > math.sqrt(hF.Integral()) : eF = math.sqrt(hF.Integral())
        rootfile.Close()

        effis['dataAltBkg'] = computeEffi(nP,nF,eP,eF)
    else:
        effis['dataAltBkg'] = [-1,-1]
        effis['canv_dataAltBkg'] = None

    return effis



def getAllScales( info, bindef, refReplica ):
    scales = {}

    for key,rfile in info.iteritems():
        if not info[key] is None and os.path.isfile(rfile) :
            rootfile = ROOT.TFile( rfile, 'read' )
            replica = int(rfile.split('_Stat')[-1].split('.root')[0]) if 'dataReplica' in key else refReplica
            fitresP = rootfile.Get( '%s_resP_Stat%d' % (bindef['name'],replica)  )
            
            fitMean = fitresP.floatParsFinal().find('meanP')
            v = fitMean.getVal()
            e = fitMean.getError()
            rootfile.Close()
        
            scales[key] = [v,e]
        else:
            scales[key] = [-999,-999]
    return scales


