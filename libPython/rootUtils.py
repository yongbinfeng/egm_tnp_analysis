#!/usr/bin/env python3 

import ROOT, math, array, ctypes, copy
import numpy as np
import os.path
from . import fitUtils
from .plotUtils import safeOpenFile, safeGetObject
import functools

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
    infile = safeOpenFile(p, mode="read")
    #print(infile.ls())
    h_tmp_pass = safeGetObject(infile, f"pass_{sample.getName()}", detach=False)
    h_tmp_fail = safeGetObject(infile, f"fail_{sample.getName()}", detach=False)
    altPass = "pass_" + sample.getName() + "_alt"
    keyNames = [k.GetName() for k in infile.GetListOfKeys()]
    h_tmp_pass_alt = None # will be needed for tracking when using all probes to form failing probe MC template to fit data, in this case the alternate version of passing probe made with  standalone variables is needed 
    if altPass in keyNames:
        h_tmp_pass_alt = safeGetObject(infile, altPass, detach=False)
        
    outfile = safeOpenFile(sample.getOutputPath(), mode="recreate")
    
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


def histPlotter(rootfile, tnpBin, plotDir, replica=-1, verbosePlotting=True ):

    if verbosePlotting:
        binName = tnpBin['name']
        if replica<0:
            c = safeGetObject(rootfile, f"{binName}_Canv", detach=False)
            c.SaveAs(f"{plotDir}/{binName}.png")   
            # c.Print( '%s/%s.pdf' % (plotDir,tnpBin['name']))  ## this doesn't work properly, it only saves the pad with the text
        else:
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
    #print("inside getAllEffi")
    binName = bindef["name"]
    for key in info.keys():
        value = info[key]
        if value is None or not os.path.isfile(value):
            # print(f"{key} -> {value}: returning empty output")
            effis[key] = [-1,-1]
            effis[f"canv_{key}"] = None
        else:
            rootfile = safeOpenFile(value, mode="read")
            # some customization
            if key == "mcNominal":
                # just counting using integral
                hP = safeGetObject(rootfile, f"{binName}_Pass", detach=False)
                hF = safeGetObject(rootfile, f"{binName}_Fail", detach=False)
                bin1 = 1
                bin2 = hP.GetXaxis().GetNbins()
                ePc = ctypes.c_double(-1.0)
                eFc = ctypes.c_double(-1.0)
                nP = hP.IntegralAndError(bin1,bin2,ePc)
                nF = hF.IntegralAndError(bin1,bin2,eFc)
                eP = float(ePc.value)
                eF = float(eFc.value)
            else:
                #print(key)
                fitresP = safeGetObject(rootfile, f"{binName}_resP", detach=False)
                fitresF = safeGetObject(rootfile, f"{binName}_resF", detach=False)
                canv = safeGetObject(rootfile, f"{binName}_Canv", detach=False)
                #effis[f"canv_{key}"] = copy.deepcopy(canv.Clone(f"{key}_{binName}"))
                effis[f"canv_{key}"] = copy.deepcopy(canv)
                fitP = fitresP.floatParsFinal().find("nSigP")
                fitF = fitresF.floatParsFinal().find("nSigF")
                nP = fitP.getVal()
                nF = fitF.getVal()
                eP = fitP.getError()
                eF = fitF.getError()
                if key in ["dataNominal", "dataAltSig", "dataAltBkg"]:
                    # if the fit is reliable the uncertainty in nSignal for data must be at least equal to sqrt(nSignal)
                    # it can be larger because of the presence of parameters and backgrounds
                    eP = max(math.sqrt(nP), eP)
                    eF = max(math.sqrt(nF), eF)
            # done with all cases
            effis[key] = computeEffi(nP,nF,eP,eF) +[nP, nF, eP, eF]
            #rootfile.Close() # try also not closing the file
    #print("Outside getAllEffi")
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


