from libPython.tnpClassUtils import tnpSample

### qll stat
#inputDir = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/tnp/2021-01-18/'
inputDir = '/data/shared/tnp/2021-02-10/'

wmass_selection = {
    # MiniAOD TnP for IDs scale factors
    'mu_DY_preVFP'        : tnpSample('mu_DY_preVFP'    , inputDir +'/dy_preVFP.root'         , isMC = True       , nEvts = -1 ) ,
    'mu_DY_postVFP'       : tnpSample('mu_DY_postVFP'   , inputDir +'/dy_postVFP.root'        , isMC = True       , nEvts = -1 ) ,

    'mu_DY_all'           : tnpSample('mu_DY_all'       , inputDir +'/dy_all.root'            , isMC = True       , nEvts = -1 ) ,

    'mu_RunB2016'         : tnpSample('mu_RunB'         , inputDir +'/data_runB.root'         , lumi =  5.75049 ) ,
    'mu_RunC2016'         : tnpSample('mu_RunC'         , inputDir +'/data_runC.root'         , lumi =  2.57649 ) ,
    'mu_RunD2016'         : tnpSample('mu_RunD'         , inputDir +'/data_runD.root'         , lumi =  4.24229 ) ,
    'mu_RunE2016'         : tnpSample('mu_RunE'         , inputDir +'/data_runE.root'         , lumi =  4.02523 ) ,
    'mu_RunF2016_preVFP'  : tnpSample('mu_RunF'         , inputDir +'/data_runF.root'         , lumi =  2.69053 ) ,
    'mu_RunF2016_postVFP' : tnpSample('mu_RunF_postVFP' , inputDir +'/data_runF_postVFP.root' , lumi =  0.41398 ) ,
    'mu_RunG2016'         : tnpSample('mu_RunG'         , inputDir +'/data_runG.root'         , lumi =  7.57582 ) ,
    'mu_RunH2016'         : tnpSample('mu_RunH'         , inputDir +'/data_runH.root'         , lumi =  8.65063 ) ,

    'mu_RunBtoF'          : tnpSample('mu_RunBtoF'      , inputDir +'/data_BtoF.root'         , lumi = 19.28503 ) ,
    'mu_RunGtoH'          : tnpSample('mu_RunGtoH'      , inputDir +'/data_GtoH.root'         , lumi = 16.64043 ) ,
    'mu_RunBtoH'          : tnpSample('mu_RunBtoH'      , inputDir +'/data_BtoH.root'         , lumi = 35.92546 ) ,

}
