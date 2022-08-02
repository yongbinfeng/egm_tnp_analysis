from libPython.tnpClassUtils import tnpSample

### qll stat
#inputDir = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/tnp/2021-01-18/'
#inputDir = '/data/shared/tnp/2021-02-10/'
#inputDir = '/data/shared/tnp/2021-03-02/'
#inputDir = '/data/shared/tnp/2021-06-07/'
#inputDir = '/data/shared/tnp/2021-07-12-newStyle/'
#inputDir = '/data/shared/tnp/2021-09-20/'
#inputDir = '/data/shared/tnp/2021-10-07-tighterGenMatchDR/'
## on the old machine inputDir = '/data/shared/tnp/2021-10-29-trackInfo/'
inputDir = '/home/users/rajarshi/Steve/Steve_Marc_Raj/' ## on the new machine


wmass_selection = {
    # MiniAOD TnP for IDs scale factors
    'mu_DY_preVFP'        : tnpSample('mu_DY_preVFP'    , inputDir +'/dy_preVFP.root'         , isMC = True       , nEvts = -1 ) ,
    'mu_DY_postVFP'       : tnpSample('mu_DY_postVFP'   , inputDir +'/Steve_MC_Isolation_20_07_2022.root'        , isMC = True       , nEvts = -1 ) ,

    'mu_DY_all'           : tnpSample('mu_DY_all'       , inputDir +'/dy_all.root'            , isMC = True       , nEvts = -1 ) ,

    'mu_RunB'         : tnpSample('mu_RunB'         , inputDir +'/data_runB.root'         , lumi =  5.75049 ) ,
    'mu_RunC'         : tnpSample('mu_RunC'         , inputDir +'/data_runC.root'         , lumi =  2.57649 ) ,
    'mu_RunD'         : tnpSample('mu_RunD'         , inputDir +'/data_runD.root'         , lumi =  4.24229 ) ,
    'mu_RunE'         : tnpSample('mu_RunE'         , inputDir +'/data_runE.root'         , lumi =  4.02523 ) ,
    'mu_RunF_preVFP'  : tnpSample('mu_RunF_preVFP'  , inputDir +'/data_runF_preVFP.root'  , lumi =  2.69053 ) ,
    'mu_RunF_postVFP' : tnpSample('mu_RunF_postVFP' , inputDir +'/data_runF_postVFP.root' , lumi =  0.41398 ) ,
    'mu_RunG'         : tnpSample('mu_RunG'         , inputDir +'/data_runG.root'         , lumi =  7.57582 ) ,
    'mu_RunH'         : tnpSample('mu_RunH'         , inputDir +'/data_runH.root'         , lumi =  8.65063 ) ,

    'mu_RunBtoF'          : tnpSample('mu_RunBtoF'      , inputDir +'/data_BtoF.root'         , lumi = 19.28503 ) ,
    'mu_RunGtoH'          : tnpSample('mu_RunGtoH'      , inputDir +'/Steve_Data_Isolation_20_07_2022.root'         , lumi = 16.64043 ) ,
    'mu_RunBtoH'          : tnpSample('mu_RunBtoH'      , inputDir +'/data_BtoH.root'         , lumi = 35.92546 ) ,

}
