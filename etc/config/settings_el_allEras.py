#############################################################
########## General settings
#############################################################
# flag to be Tested
cutpassTrigEl = '( probe_eleTrgPt  > -1 )'

# flag to be Tested
flags = {
    'triggerEl'  : cutpassTrigEl,
    }
baseOutDir = 'results/elFullData/'

#############################################################
########## samples definition  - preparing the samples
#############################################################
### samples are defined in etc/inputs/tnpSampleDef.py
### not: you can setup another sampleDef File in inputs
import etc.inputs.tnpSampleDef as tnpSamples
tnpTreeDir = 'IDIsoToHLT'

samplesDef = {
    'data'   : tnpSamples.wmass_80X['el_Run2016B'].clone(),
    'mcNom'  : tnpSamples.wmass_80X['el_DY'].clone(),
    'mcAlt'  : None, #' #tnpSamples.ICHEP2016['mc_DY_amcatnlo_ele'].clone(),
    'tagSel' : None, #tnpSamples.ICHEP2016['mc_DY_madgraph_ele'].clone(),
}
## can add data sample easily
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016C'] )
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016D'] )
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016E'] )
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016F'] )
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016G'] )
samplesDef['data'].add_sample( tnpSamples.wmass_80X['el_Run2016H'] )

## some sample-based cuts... general cuts defined here after
## require mcTruth on MC DY samples and additional cuts
## all the samples MUST have different names (i.e. sample.name must be different for all)
## if you need to use 2 times the same sample, then rename the second one
#samplesDef['data'  ].set_cut('run >= 273726')

## marc if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_mcTruth()
## marc if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_mcTruth()
## marc if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_mcTruth()
## marc if not samplesDef['tagSel'] is None:
## marc     samplesDef['tagSel'].rename('mcAltSel_DY_madgraph_ele')
## marc     samplesDef['tagSel'].set_cut('tag_Ele_pt > 33  && tag_Ele_nonTrigMVA > 0.90')

## set MC weight, simple way (use tree weight) 
weightName = '1.'
## marc if not samplesDef['mcNom' ] is None: samplesDef['mcNom' ].set_weight(weightName)
## marc if not samplesDef['mcAlt' ] is None: samplesDef['mcAlt' ].set_weight(weightName)
## marc if not samplesDef['tagSel'] is None: samplesDef['tagSel'].set_weight(weightName)

#############################################################
########## bining definition  [can be nD bining]
#############################################################
binning_eta = [-2.5,-2.3,-2.1,-1.9,-1.7,-1.566,-1.4442,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]
binning_pt  = [30, 31.5, 33, 36, 39, 42, 45]

biningDef = [
   { 'var' : 'probe_lep_eta', 'type': 'float', 'bins': binning_eta },
   { 'var' : 'probe_lep_pt' , 'type': 'float', 'bins': binning_pt  },
]

#############################################################
########## Cuts definition for all samples
#############################################################
### cut
cutBase   = 'tag_lep_pt > 30 ' ## cuts already applied in the tree

# can add addtionnal cuts for some bins (first check bin number using tnpEGM --checkBins)
## additionalCuts = { 
##     0 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     1 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     2 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     3 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     4 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     5 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     6 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     7 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     8 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45',
##     9 : 'tag_Ele_trigMVA > 0.92 && sqrt( 2*event_met_pfmet*tag_Ele_pt*(1-cos(event_met_pfphi-tag_Ele_phi))) < 45'
## }

#### or remove any additional cut (default)
additionalCuts = None

#############################################################
########## fitting params to tune fit by hand if necessary
#############################################################
tnpParNomFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
    "acmsP[60.,50.,80.]","betaP[0.05,0.01,0.08]","gammaP[0.1, 0, 1]","peakP[90.0]",
    "acmsF[60.,50.,80.]","betaF[0.05,0.01,0.08]","gammaF[0.1, 0, 1]","peakF[90.0]",
    ]

tnpParAltSigFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[1,0.7,6.0]","alphaP[2.0,1.2,3.5]" ,'nP[3,-5,5]',"sigmaP_2[1.5,0.5,6.0]","sosP[1,0.5,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[2,0.7,15.0]","alphaF[2.0,1.2,3.5]",'nF[3,-5,5]',"sigmaF_2[2.0,0.5,6.0]","sosF[1,0.5,5.0]",
    "acmsP[60.,50.,75.]","betaP[0.04,0.01,0.06]","gammaP[0.1, 0.005, 1]","peakP[90.0]",
    "acmsF[60.,50.,75.]","betaF[0.04,0.01,0.06]","gammaF[0.1, 0.005, 1]","peakF[90.0]",
    ]
     
tnpParAltBkgFit = [
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
    "alphaP[0.,-5.,5.]",
    "alphaF[0.,-5.,5.]",
    ]
        
