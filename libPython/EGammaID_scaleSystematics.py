#!/usr/bin/env python

import sys,os
from array import array
from math import sqrt
import numpy as np
import ROOT as rt
import CMS_lumi, tdrstyle

from EGammaID_scaleFactors import isFloat
from efficiencyUtils import efficiency
from efficiencyUtils import efficiencyList
import efficiencyUtils as effUtil

tdrstyle.setTDRStyle()

def doEGM_Scales(filein, lumi, axis = ['p_{T} (GeV)','#eta'] ):
    print(" Opening file: {f} (plot lumi: {l:.1f})".format(f=filein,l=lumi ))
    CMS_lumi.lumi_13TeV = "%3.1f fb^{-1}" % lumi 

    nameOutBase = filein 
    if not os.path.exists( filein ) :
        print('file {f} does not exist'.format(f=filein))
        sys.exit(1)


    fileWithScales = open(filein, 'r')

    systs = {}
    etaEdges = []; ptEdges = []    
    for line in fileWithScales :
        modifiedLine = line.lstrip(' ').rstrip(' ').rstrip('\n')
        numbers = modifiedLine.split('\t')        

        if len(numbers) > 0 and isFloat(numbers[0]):
            nominal = float(numbers[4])
            scales = numbers[5:]
            etaKey = ( float(numbers[0]), float(numbers[1]) )
            ptKey  = ( float(numbers[2]), float(numbers[3]) ) 
            if float(numbers[0]) not in etaEdges: etaEdges.append(float(numbers[0]))
            if float(numbers[2]) not in ptEdges:  ptEdges.append(float(numbers[2]))
            
            mysysts = [(float(scales[i])-nominal)/91.18/sqrt(2.) for i in range(len(scales))]
            systs[(etaKey,ptKey)] = mysysts

    etaEdges.append(float(numbers[1]))
    ptEdges.append(float(numbers[3]))

    fileWithScales.close()
    
    systHistos2D = []
    print("etaEdges = ",etaEdges)
    print("ptEdges = ",ptEdges)
    
    systHisto = rt.TH2F('systHisto','',len(etaEdges)-1,array('f',etaEdges),len(ptEdges)-1,array('f',ptEdges))
    systHisto.GetYaxis().SetTitle(axis[0])
    systHisto.GetXaxis().SetTitle(axis[1])
    for isyst in range(len(scales)):
        if isyst==0: # this is very large above 48 GeV
            systHisto.GetZaxis().SetRangeUser(-2e-2,2e-2)
        else: # these are statistical replicas
            systHisto.GetZaxis().SetRangeUser(-2e-3,2e-3)
        name = ('syst_%d' % isyst) if isyst<2 else ('stat_%d' % (isyst-2))
        thish = systHisto.Clone(name)
        systHistos2D.append(thish)

    for (etaKey,ptKey),systValues in systs.iteritems():
        eta = np.mean(etaKey)
        pt = np.mean(ptKey)
        etabin = systHisto.GetXaxis().FindBin(eta)
        ptbin  = systHisto.GetYaxis().FindBin(pt)
        for i in range(len(systValues)):
            systHistos2D[i].SetBinContent(etabin,ptbin,systValues[i])

    pdfout = nameOutBase + '_egammaPlots.pdf'
    cDummy = rt.TCanvas('','',900,600)
    cDummy.Print( pdfout + "[" )

    rt.gStyle.SetPalette(51)
    rt.gStyle.SetPaintTextFormat('1.4f');
    rt.TColor.CreateGradientColorTable(3,
                                      array ("d", [0.00, 0.50, 1.00]),
                                      array ("d", [0.00, 1.00, 1.00]),
                                      array ("d", [0.34, 1.00, 0.65]),
                                      array ("d", [0.82, 1.00, 0.00]),
                                      255,  0.95)

    c2D = rt.TCanvas('canScaleFactor','canScaleFactor',900,600)
    c2D.SetRightMargin(0.15)
    c2D.SetLeftMargin( 0.15)
    c2D.SetTopMargin(  0.10)
    for ih in systHistos2D:
        c2D.SetTitle(ih.GetName())
        ih.DrawCopy("colz0 TEXT45")
        c2D.Print( pdfout )

    cDummy.Print( pdfout + "]" )

    rootout = rt.TFile(nameOutBase + '_EGM2D.root','recreate')
    rootout.cd()
    for ih in systHistos2D:
        ih.Write()
    rootout.Close()

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description='tnp EGM scale factors')
    parser.add_argument('--lumi'  , type = float, default = -1, help = 'Lumi (just for plotting purpose)')
    parser.add_argument('--txtFile' , default = None, help = 'EGM formatted txt file')
    args = parser.parse_args()

    if args.txtFile is None:
        print(' - Needs EGM txt file as input')
        sys.exit(1)
    

    CMS_lumi.lumi_13TeV = "5.5 fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.lumi_sqrtS = "13 TeV"
    
    axis = ['p_{T} (GeV)','#eta']

    doEGM_Scales(args.txtFile, args.lumi,axis)
