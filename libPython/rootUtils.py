import ROOT, math, array
import numpy as np
from fitUtils import *
import functools
#from fitSimultaneousUtils import *

def removeNegativeBins(h):
    for i in xrange(h.GetNbinsX()):
	if (h.GetBinContent(i) < 0):
            h.SetBinContent(i, 0)


def makePassFailHistograms( sample, flag, bins, bindef, commonCuts, var ):
    ## open rootfile
    tree = rt.TChain(sample.tree)
    for p in sample.path:
        print ' adding rootfile: ', p
        tree.Add(p)
    
    if not sample.puTree is None:
        print ' - Adding weight tree: %s from file %s ' % (sample.weight.split('.')[0], sample.puTree)
        tree.AddFriend(sample.weight.split('.')[0],sample.puTree)


    ## open outputFile
    hPass = []
    hFail = []
    lCuts = []

    notflag = '!({f})'.format(f=flag)

    ## now slowly constructing the 3d histogram. some things here are hardcoded...

    probe_binning_eta, probe_binning_pt = (bindef[0]['bins'], bindef[1]['bins']) if 'eta' in bindef[0]['var'] else (bindef[1]['bins'], bindef[0]['bins'])
    probe_var_eta, probe_var_pt         = (bindef[0]['var'] , bindef[1]['var'] ) if 'eta' in bindef[0]['var'] else (bindef[1]['var'] , bindef[0]['var'] )

    probe_binning_pt  = array.array('d', probe_binning_pt)
    probe_binning_eta = array.array('d', probe_binning_eta)

    binning_mass = array.array('d', [var['min'] + i*(var['max']-var['min'])/var['nbins'] for i in range(var['nbins']+1)])

    ## now have the full 3d binning


    cuts = commonCuts
    if sample.mcTruth :
        cuts = ' && mcTrue==1'
    if not sample.cut is None :
        cuts = ' && {sc}'.format(sc=sample.cut)

    ## get the weight to fill the histogram with

    if not sample.weight is None:
        tmp_weight = sample.weight
        print 'found sample weight of', sample.weight
    else:
        tmp_weight = 1.
        print 'did not find sample weight, setting to 1.'

    if sample.maxWeight < 999:
        tmp_weight = '({w} < {mw} ? {w} : 1.0 )' .format(w=sample.weight,mw=sample.maxWeight)

    cutPass = '( {c} && {f} )' .format(c=cuts, f=   flag)
    cutFail = '( {c} && {f} )' .format(c=cuts, f=notflag)

    ##fill the passing and failing 3D histogram
    h_tmp_pass =  ROOT.TH3D('pass_'+sample.name, 'pass_'+sample.name, var['nbins']            , binning_mass, 
                                                                      len(probe_binning_pt) -1, probe_binning_pt ,
                                                                      len(probe_binning_eta)-1, probe_binning_eta)
    h_tmp_fail =  ROOT.TH3D('fail_'+sample.name, 'fail_'+sample.name, var['nbins']            , binning_mass, 
                                                                      len(probe_binning_pt) -1, probe_binning_pt ,
                                                                      len(probe_binning_eta)-1, probe_binning_eta)
    h_tmp_pass.Sumw2()
    h_tmp_fail.Sumw2()

    print 'will fill 3d histogram with:' 
    print '   cut:', cuts
    print '   using weight', tmp_weight
    print '   passing flag', flag

    print 'now filling passing histogram for sample', sample.name
    drawret = tree.Draw('{z}:{y}:{x}>>{h}'.format(z=probe_var_eta,y=probe_var_pt,x=var['name'],h=h_tmp_pass.GetName()), '({c})*({w})'.format(c=cutPass,w=tmp_weight))
    if drawret == -1:
        print 'some error occured, please check!'
        exit(0)

    print 'now filling failing histogram for sample', sample.name
    drawret = tree.Draw('{z}:{y}:{x}>>{h}'.format(z=probe_var_eta,y=probe_var_pt,x=var['name'],h=h_tmp_fail.GetName()), '({c})*({w})'.format(c=cutFail,w=tmp_weight))
    if drawret == -1:
        print 'some error occured, please check!'
        exit(0)

    outfile = ROOT.TFile(sample.histFile,'recreate')

    for ii,ib in enumerate(bins):
        h_name = ib['name' ]
        h_title= ib['title']

        tmp_valpt  = ib['vars'][probe_var_pt ]['min']
        tmp_valeta = ib['vars'][probe_var_eta]['min']

        ibin_pt  = h_tmp_pass.GetYaxis().FindBin(tmp_valpt )
        ibin_eta = h_tmp_pass.GetZaxis().FindBin(tmp_valeta)
        #print 'i am at ibin_pt  {ipt}  and pt  {pt:.1f} '.format(ipt = ibin_pt , pt = tmp_valpt )
        #print 'i am at ibin_eta {ieta} and eta {eta:.1f}'.format(ieta= ibin_eta, eta= tmp_valeta)

        h_pass = h_tmp_pass.ProjectionX(h_name+'_Pass', ibin_pt, ibin_pt, ibin_eta, ibin_eta)
        h_fail = h_tmp_fail.ProjectionX(h_name+'_Fail', ibin_pt, ibin_pt, ibin_eta, ibin_eta)
        h_pass .SetTitle(h_title+' passing')
        h_fail .SetTitle(h_title+' failing')

        removeNegativeBins(h_pass)
        removeNegativeBins(h_fail)

        h_pass   .Write(h_pass   .GetName())
        h_fail   .Write(h_fail   .GetName())

        bin1 = 1
        bin2 = h_pass.GetXaxis().GetNbins()
        epass = ROOT.Double(-1.0)
        efail = ROOT.Double(-1.0)
        passI = h_pass.IntegralAndError(bin1,bin2,epass)
        failI = h_fail.IntegralAndError(bin1,bin2,efail)
        eff   = 0
        e_eff = 0
        if passI > 0 :
            itot  = (passI+failI)
            eff   = passI / (passI+failI)
            e_eff = math.sqrt(passI*passI*efail*efail + failI*failI*epass*epass) / (itot*itot)
        print ib['cut']
        print '    ==> pass: %.1f +/- %.1f ; fail : %.1f +/- %.1f : eff: %1.3f +/- %1.3f' % (passI,epass,failI,efail,eff,e_eff)

    outfile.Close()


def makeBootstrapHistograms( sample, flag, bindef, var, resample ):
    ## open rootfile
    tree = ROOT.TChain(sample.tree)
    for p in sample.path:
        print ' adding rootfile: ', p
        tree.Add(p)
    
    if not sample.puTree is None:
        print ' - Adding weight tree: %s from file %s ' % (sample.weight.split('.')[0], sample.puTree)
        tree.AddFriend(sample.weight.split('.')[0],sample.puTree)

    ## open outputFile
    outfilename = getattr(sample,'histFile{ir}'.format(ir=resample))
    print "histograms output file = ",outfilename
    outfile = ROOT.TFile(outfilename,'recreate')

    seed = 123456789+resample
    np.random.seed(seed)
    for ib in range(len(bindef['bins'])):

        ## select the events passing cuts
        cuts = bindef['bins'][ib]['cut']
        if sample.mcTruth :
            cuts = '%s && mcTrue==1' % cuts
        if not sample.cut is None :
            cuts = '%s && %s' % (cuts,sample.cut)
        
        notflag = '!(%s)' % flag
        
        if sample.isMC and not sample.weight is None:
            cutPass = '( %s && %s ) * %s ' % (cuts,    flag, sample.weight)
            if sample.maxWeight < 999:
                cutPass = '( %s && %s ) * (%s < %f ? %s : 1.0 )' % (cuts,    flag, sample.weight,sample.maxWeight,sample.weight)
        else:
            cutPass = '( %s && %s )' % (cuts,    flag)
     
        tree.Draw('>>elist',cutPass)
        elist = ROOT.gDirectory.Get('elist')

        print "Tot events = ",tree.GetEntries()," selected by the cut ",cutPass," = ",elist.GetN()

        print "Resampling # ",resample
        ## get the list of resampled events
        entriesList = range(elist.GetN())
        resamples = np.random.choice(entriesList, size=len(entriesList))
        
        hPass = ROOT.TH1D('{name}_Stat{i}'.format(name=bindef['bins'][ib]['name'],i=resample),bindef['bins'][ib]['title'],var['nbins'],var['min'],var['max'])
        hPass.Sumw2()
    
        ## fill the histograms
        print "Now looping on the resampled dataset to fill ",hPass.GetName()
        tree.SetBranchStatus("*",0)
        tree.SetBranchStatus(var['name'],1)
        for ie,entry in enumerate(resamples):
            if ie%1000==0: print "Processing selected event ",ie," / ",len(resamples)
            tev = elist.GetEntry(entry)
            tree.GetEntry(tev)
            hPass.Fill(getattr(tree,var['name']))
        tree.SetBranchStatus("*",1)

        removeNegativeBins(hPass)
     
        hPass.Write(hPass.GetName())
     
        bin1 = 1
        bin2 = hPass.GetXaxis().GetNbins()
        epass = ROOT.Double(-1.0)
        passI = hPass.IntegralAndError(bin1,bin2,epass)
        print cuts
        print '    ==> Resample ',resample,' pass: %.1f +/- %.1f ' % (passI,epass)
    outfile.Close()



def histPlotter( filename, tnpBin, plotDir, replica=-1 ):
    print 'opening ', filename
    rootfile = ROOT.TFile(filename,"read")

    if replica<0:
        print '  get canvas: ' , '%s_Canv' % tnpBin['name']
        c = rootfile.Get( '%s_Canv' % tnpBin['name'] )
        c.Print( '%s/%s.png' % (plotDir,tnpBin['name']))
        c.Print( '%s/%s.pdf' % (plotDir,tnpBin['name']))
    else:
        print '  get canvas: ' , '%s_Canv_Stat%d' % (tnpBin['name'],replica)
        c = rootfile.Get( '%s_Canv_Stat%d' % (tnpBin['name'],replica) )
        c.Print( '%s/%s_Stat%d.png' % (plotDir,tnpBin['name'],replica))
        c.Print( '%s/%s_Stat%d.pdf' % (plotDir,tnpBin['name'],replica))
        


def computeEffi( n1,n2,e1,e2):
    effout = []
    if (n1+n2):
        eff   = n1/(n1+n2)
        e_eff = 1/(n1+n2)*math.sqrt(e1*e1*n2*n2+e2*e2*n1*n1)/(n1+n2)
        if e_eff < 0.001 : e_eff = 0.001
    else:
        eff, e_eff = 1.1, 0.01

    effout.append(eff)
    effout.append(e_eff)
    
    return effout


import os.path
def getAllEffi( info, bindef ):
    effis = {}
    if not info['mcNominal'] is None and os.path.isfile(info['mcNominal']):
        rootfile = ROOT.TFile( info['mcNominal'], 'read' )
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])
        bin1 = 1
        bin2 = hP.GetXaxis().GetNbins()
        eP = ROOT.Double(-1.0)
        eF = ROOT.Double(-1.0)
        nP = hP.IntegralAndError(bin1,bin2,eP)
        nF = hF.IntegralAndError(bin1,bin2,eF)

        effis['mcNominal'] = computeEffi(nP,nF,eP,eF)
        rootfile.Close()
    else: effis['mcNominal'] = [-1,-1]

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
        hP = rootfile.Get('%s_Pass'%bindef['name'])
        hF = rootfile.Get('%s_Fail'%bindef['name'])
        bin1 = 1
        bin2 = hP.GetXaxis().GetNbins()
        eP = ROOT.Double(-1.0)
        eF = ROOT.Double(-1.0)
        nP = hP.IntegralAndError(bin1,bin2,eP)
        nF = hF.IntegralAndError(bin1,bin2,eF)

        effis['mcAlt'] = computeEffi(nP,nF,eP,eF)
        rootfile.Close()
    else: effis['mcAlt'] = [-1,-1]

    if not info['dataNominal'] is None and os.path.isfile(info['dataNominal']) :
        rootfile = ROOT.TFile( info['dataNominal'], 'read' )
        from ROOT import RooFit,RooFitResult
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )

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

        effis['dataNominal'] = computeEffi(nP,nF,eP,eF)
    else:
        effis['dataNominal'] = [-1,-1]
    if not info['dataAltSig'] is None and os.path.isfile(info['dataAltSig']) :
        rootfile = ROOT.TFile( info['dataAltSig'], 'read' )
        from ROOT import RooFit,RooFitResult
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )

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

        effis['dataAltSig'] = computeEffi(nP,nF,eP,eF)

    else:
        effis['dataAltSig'] = [-1,-1]

    if not info['dataAltBkg'] is None and os.path.isfile(info['dataAltBkg']):
        rootfile = ROOT.TFile( info['dataAltBkg'], 'read' )
        from ROOT import RooFit,RooFitResult
        fitresP = rootfile.Get( '%s_resP' % bindef['name']  )
        fitresF = rootfile.Get( '%s_resF' % bindef['name'] )

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
    return effis



def getAllScales( info, bindef, refReplica ):
    scales = {}

    for key,rfile in info.iteritems():
        if not info[key] is None and os.path.isfile(rfile) :
            rootfile = ROOT.TFile( rfile, 'read' )
            replica = int(rfile.split('_Stat')[-1].split('.root')[0]) if 'dataReplica' in key else refReplica
            from ROOT import RooFit,RooFitResult
            fitresP = rootfile.Get( '%s_resP_Stat%d' % (bindef['name'],replica)  )
            
            fitMean = fitresP.floatParsFinal().find('meanP')
            v = fitMean.getVal()
            e = fitMean.getError()
            rootfile.Close()
        
            scales[key] = [v,e]
        else:
            scales[key] = [-999,-999]
    return scales


