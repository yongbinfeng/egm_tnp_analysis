from libPython.tnpClassUtils import tnpSample

### qll stat
eosDirTrigger = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/tnp/2020-10-28-firstTest/RunG/'

wmass_80X = {
    # MiniAOD TnP for IDs scale factors
    'mu_DY' : tnpSample('mu_DY', eosDirTrigger+'/triggerTnP_DYJetsToLL_M50_all_mu.root', isMC = True, nEvts = -1 ),

    'mu_Run2016_all' : tnpSample('mu_Run2016_all' , eosDirTrigger+'/triggerTnP_SingleMuon_Run2016_all.root' , lumi = 35.93 ),
    'mu_Run2016_all_othertrees' : tnpSample('mu_Run2016_all_othertrees' , '/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_recoLeptons_and_triggerMatch_latest/triggerTnP_SingleMuon_Run2016_all.root' , lumi = 35.93 ),

    'mu_DY_noMu50'         : tnpSample('mu_DY_noMu50'          , '/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_TNP_2l_DY_loosestLep_triggerMatch/triggerTnP_DY_mu.root', isMC = True, nEvts = -1 ),

    'mu_Run2016_noMu50_all': tnpSample('mu_Run2016_noMu50_all' , '/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_TNP_2l_SingleMu_loosestLep_triggerMatch/triggerTnP_SingleMuon_Run2016_all_mu.root' , lumi = 35.93 ),
    'mu_Run2016_noMu50_G': tnpSample('mu_Run2016_noMu50_G' , '/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_TNP_2l_SingleMu_loosestLep_triggerMatch/triggerTnP_SingleMuon_Run2016G_mu.root' , lumi = 7.58 ),
    'mu_Run2016_noMu50_H': tnpSample('mu_Run2016_noMu50_H' , '/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_TNP_2l_SingleMu_loosestLep_triggerMatch/triggerTnP_SingleMuon_Run2016H_mu.root' , lumi = 8.65 ),

## 
##     # MiniAOD TnP for IDs scale factors
     'el_DY'          : tnpSample('el_DY'         , eosDirTrigger+'/triggerTnP_DYJetsToLL_M50_all_el.root'     , isMC = True, nEvts = -1 ),
## 
     'el_Run2016_all' : tnpSample('el_Run2016_all', eosDirTrigger+'/triggerTnP_SingleElectron_Run2016_all.root', lumi = 35.93 ),
}

### qll stat
eosDirSelection = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/tnp/2020-11-11-withTrigMatch/'

wmass_selection = {
    # MiniAOD TnP for IDs scale factors
    'mu_DY_noAPV'   : tnpSample('mu_DY_noAPV' , eosDirSelection +'/DY_noAPV.root'  , isMC = True, nEvts = -1 ),
    'mu_DY_APV'   : tnpSample('mu_DY_APV' , eosDirSelection +'/DYAPV.root'  , isMC = True, nEvts = -1 ),

    'mu_RunG2016'   : tnpSample('mu_RunG_all' , eosDirSelection +'/data_RunG.root'  , lumi =  7.57582 ),
    'mu_RunH2016'   : tnpSample('mu_RunH_all' , eosDirSelection +'/data_RunH.root'  , lumi =  8.65063 ),
    'mu_RunGH2016'  : tnpSample('mu_RunGH_all', eosDirSelection +'/data_RunGH.root' , lumi = 16.22645 ),

}
