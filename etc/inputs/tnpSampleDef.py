from libPython.tnpClassUtils import tnpSample

### qll stat
#eosDirElData = '/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_ELECTRONS/'
eosDirElData = '/afs/cern.ch/work/e/emanuele/wmass/tnp/CMSSW_10_2_0/src/egm_tnp_analysis/localtrees/'
eosDirMuData = '/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_MUONS/'
eosDirMC     = '/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_MC/'

wmass_80X = {
    # MiniAOD TnP for IDs scale factors
    'mu_DY' : tnpSample('mu_DY', eosDirMC+'/triggerTnP_DYJetsToLL_M50_all_mu.root', isMC = True, nEvts = -1 ),

    'mu_Run2016B' : tnpSample('mu_Run2016B' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016B_mu.root' , lumi = 5.767 ),
    'mu_Run2016C' : tnpSample('mu_Run2016C' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016C_mu.root' , lumi = 2.646 ),
    'mu_Run2016D' : tnpSample('mu_Run2016D' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016D_mu.root' , lumi = 4.353 ),
    'mu_Run2016E' : tnpSample('mu_Run2016E' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016E_mu.root' , lumi = 3.985 ),
    'mu_Run2016F' : tnpSample('mu_Run2016F' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016F_mu.root' , lumi = 3.160 ),
    'mu_Run2016G' : tnpSample('mu_Run2016G' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016G_mu.root' , lumi = 7.539 ),
    'mu_Run2016H' : tnpSample('mu_Run2016H' , eosDirMuData+'/triggerTnP_SingleMuon_Run2016H_mu.root' , lumi = 8.762 ),

    # MiniAOD TnP for IDs scale factors
    'el_DY' : tnpSample('mu_DY', eosDirMC+'/triggerTnP_DYJetsToLL_M50_all_el.root', isMC = True, nEvts = -1 ),

    'el_Run2016' : tnpSample('el_Run2016' , eosDirElData+'/triggerTnP_SingleElectron_Run2016_el.root' , lumi = 36.212 ),

    # 'el_Run2016B' : tnpSample('el_Run2016B' , eosDirElData+'/triggerTnP_SingleElectron_Run2016B_el.root' , lumi = 5.767 ),
    # 'el_Run2016C' : tnpSample('el_Run2016C' , eosDirElData+'/triggerTnP_SingleElectron_Run2016C_el.root' , lumi = 2.646 ),
    # 'el_Run2016D' : tnpSample('el_Run2016D' , eosDirElData+'/triggerTnP_SingleElectron_Run2016D_el.root' , lumi = 4.353 ),
    # 'el_Run2016E' : tnpSample('el_Run2016E' , eosDirElData+'/triggerTnP_SingleElectron_Run2016E_el.root' , lumi = 3.985 ),
    # 'el_Run2016F' : tnpSample('el_Run2016F' , eosDirElData+'/triggerTnP_SingleElectron_Run2016F_el.root' , lumi = 3.160 ),
    # 'el_Run2016G' : tnpSample('el_Run2016G' , eosDirElData+'/triggerTnP_SingleElectron_Run2016G_el.root' , lumi = 7.539 ),
    # 'el_Run2016H' : tnpSample('el_Run2016H' , eosDirElData+'/triggerTnP_SingleElectron_Run2016H_el.root' , lumi = 8.762 ),
}
