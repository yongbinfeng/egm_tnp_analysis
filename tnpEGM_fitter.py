
### python specific import
import os
import sys
import pickle
import shutil
from multiprocessing import Pool

def joinCuts(_list):
    return '(' + ' && '.join(_list) + ')'

from optparse import OptionParser
#parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
#parser.add_option('-c', '--channel'   , dest='channel'   , type='string'      , default='mu',  help='run tnp ntuples for muons/electrons. default mu')
#parser.add_option('-i', '--indir'     , dest='indir'     , type='string'      , default=''  ,  help='directory with the trees')


parser = OptionParser(usage="%prog [options] ")
parser.add_option('--checkBins'  , action='store_true'  , help = 'check  bining definition')
parser.add_option('--createBins' , action='store_true'  , help = 'create bining definition')
parser.add_option('--createHists', action='store_true'  , help = 'create histograms')
parser.add_option('--sample'     , default='all'        , help = 'create histograms (per sample, expert only)')
parser.add_option('--altSig'     , action='store_true'  , help = 'alternate signal model fit')
parser.add_option('--altBkg'     , action='store_true'  , help = 'alternate background model fit')
parser.add_option('--doFit'      , action='store_true'  , help = 'fit sample (sample should be defined in settings.py)')
parser.add_option('--mcSig'      , action='store_true'  , help = 'fit MC nom [to init fit parama]')
parser.add_option('--doPlot'     , action='store_true'  , help = 'plotting')
parser.add_option('--sumUp'      , action='store_true'  , help = 'sum up efficiencies')
parser.add_option('--iBin'       , dest = 'binNumber'   , type = int,  default=-1, help='bin number (to refit individual bin)')
parser.add_option('--flag'       , default = None       , help ='WP to test')
parser.add_option('--era'        , type='string', default = '', help ='era to perform tnp fits for. options are: [BtoF, GtoH, BtoH]')
##parser.add_argument('settings'     , default = None       , help = 'setting file [mandatory]')


(options, args) = parser.parse_args()
#args = parser.parse_args()

##print('===> settings %s <===' % args.settings)
##importSetting = 'import %s as tnpConf' % args.settings.replace('/','.').split('.py')[0]
##print(importSetting)
##exec(importSetting)

### tnp library
import libPython.binUtils  as tnpBiner
import libPython.rootUtils as tnpRoot
import libPython.fitUtils as fitUtils
#from libPython import rootUtils as tnpRoot


if options.flag is None:
    print('[tnpEGM_fitter] flag is MANDATORY, this is the working point as defined in the settings.py')
    sys.exit(0)
    
## put most of the stuff here and make it configurable...
## ===========================================================================


## define the binning here, much easier...
binning_eta = [-2.4+0.1*i for i in range(49) ]
#binning_eta = [-2.4, -2.25, -2.10, -1.95, -1.8, -1.7, -1.:q

binning_pt  = [25., 27.5, 30., 32., 34, 36., 38., 40., 42., 44., 47., 50., 55., 65]

binningDef = [
   { 'var' : 'probe_eta', 'type': 'float', 'bins': binning_eta },
   { 'var' : 'probe_pt' , 'type': 'float', 'bins': binning_pt  },
]

weightName = 'totWeight/std::abs(totWeight)*std::max(4.f,std::abs(totWeight))'
#weightName = 'totWeight/TMath::Abs(totWeight)*TMath::Max(4.,TMath::Abs(totWeight))'

## now the flags and the cuts, pretty straight forward
cutBase     = 'tag_pt > 25. && probe_pt > 25.'

cutMinus    = 'probe_charge < 0'
cutPlus     = 'probe_charge > 0'
cutTracking = '(probe_isGlobal > 0.5 || probe_isTracker > 0.5)'
cutTrigger  = 'probe_triggerMatch > 0.5'
cutIdIp     = 'probe_mediumId > 0.5 && fabs(probe_dxy) < 0.05 && fabs(probe_dz) < 0.2'
cutIso      = 'probe_iso < 0.15'

## any additional cut to apply
additionalCuts = None

## the first item in the flags is the flag to test, the second is the base cut to be applied to all events
flags = {
    'mu_tracking_both'  : (cutTracking , joinCuts([cutBase          ]) ),
    'mu_tracking_minus' : (cutTracking , joinCuts([cutBase, cutMinus]) ),
    'mu_tracking_plus'  : (cutTracking , joinCuts([cutBase, cutPlus ]) ),

    'mu_idip_both'      : (cutIdIp     , joinCuts([cutBase,           cutTracking]) ),
    'mu_idip_minus'     : (cutIdIp     , joinCuts([cutBase, cutMinus, cutTracking]) ),
    'mu_idip_plus'      : (cutIdIp     , joinCuts([cutBase, cutPlus , cutTracking]) ),

    'mu_trigger_both'   : (cutTrigger  , joinCuts([cutBase,           cutTracking, cutIdIp]) ),
    'mu_trigger_minus'  : (cutTrigger  , joinCuts([cutBase, cutMinus, cutTracking, cutIdIp]) ),
    'mu_trigger_plus'   : (cutTrigger  , joinCuts([cutBase, cutPlus , cutTracking, cutIdIp]) ),

    'mu_iso_both'       : (cutIso      , joinCuts([cutBase,           cutTracking, cutIdIp, cutTrigger]) ),
    'mu_iso_minus'      : (cutIso      , joinCuts([cutBase, cutMinus, cutTracking, cutIdIp, cutTrigger]) ),
    'mu_iso_plus'       : (cutIso      , joinCuts([cutBase, cutPlus , cutTracking, cutIdIp, cutTrigger]) ),

    'mu_isonotrig_both' : (cutIso      , joinCuts([cutBase,           cutTracking, cutIdIp]) ),
    'mu_isonotrig_minus': (cutIso      , joinCuts([cutBase, cutMinus, cutTracking, cutIdIp]) ),
    'mu_isonotrig_plus' : (cutIso      , joinCuts([cutBase, cutPlus , cutTracking, cutIdIp]) )
}


import etc.inputs.tnpSampleDef as tnpSamples

samples_data_preVFP  = tnpSamples.wmass_selection['mu_RunBtoF'  ].clone()
samples_data_postVFP = tnpSamples.wmass_selection['mu_RunGtoH'  ].clone()
samples_data_all     = tnpSamples.wmass_selection['mu_RunBtoH'  ].clone()

samples_dy_preVFP  = tnpSamples.wmass_selection['mu_DY_preVFP' ].clone()
samples_dy_postVFP = tnpSamples.wmass_selection['mu_DY_postVFP'].clone()
samples_dy_all     = tnpSamples.wmass_selection['mu_DY_all'].clone()


tnpTreeDir = 'IDIsoToHLT'

if options.era == 'BtoF':
    samples_data = samples_data_preVFP
    samples_dy   = samples_dy_preVFP
elif options.era == 'GtoH':
    samples_data = samples_data_postVFP
    samples_dy   = samples_dy_postVFP
elif options.era == 'BtoH':
    samples_data = samples_data_all
    samples_dy   = samples_dy_all
else:
    print('you gave a wrong era name. options are \'BtoF\', \'GtoH\', and \'BtoH\'')
    sys.exit(1)

samplesDef = {
    'data'   : samples_data,
    'mcNom'  : samples_dy,
    'mcAlt'  : None, #tnpSamples.wmass_80X['mu_DY'].clone(), #tnpSamples.ICHEP2016['mc_DY_amcatnlo_ele'].clone(),
    'tagSel' : None, #tnpSamples.ICHEP2016['mc_DY_madgraph_ele'].clone(),
}

if not samplesDef['mcNom' ] == None: samplesDef['mcNom' ].set_weight(weightName)
if not samplesDef['mcAlt' ] == None: samplesDef['mcAlt' ].set_weight(weightName)
if not samplesDef['tagSel'] == None: samplesDef['tagSel'].set_weight(weightName)

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

baseOutDir = 'results/efficiencies_{era}/'.format(era=options.era)

## done making it more configurable
## ===========================================================================================

if not options.flag in flags.keys() :
    print('[tnpEGM_fitter] flag {f} not found in flags definitions'.format(f=options.flag))
    print('  --> define in settings first')
    print('  In settings I found flags: ')
    print(flags.keys())
    sys.exit(1)


outputDirectory = '%s/%s/' % (baseOutDir,options.flag)

print('===>  Output directory: ')
print(outputDirectory)


####################################################################
##### Create (check) Bins
####################################################################
if options.checkBins:
    tnpBins = tnpBiner.createBins(binningDef,flags[options.flag][1])
    tnpBiner.tuneCuts( tnpBins, additionalCuts )
    for ib in range(len(tnpBins['bins'])):
        print(tnpBins['bins'][ib]['name'])
        print('  - cut: ',tnpBins['bins'][ib]['cut'])
    sys.exit(0)
    
if options.createBins:
    if os.path.exists( outputDirectory ):
        shutil.rmtree( outputDirectory )
    os.makedirs( outputDirectory )
    tnpBins = tnpBiner.createBins(binningDef,flags[options.flag][1])
    tnpBiner.tuneCuts( tnpBins, additionalCuts )
    pickle.dump( tnpBins, open( '%s/bining.pkl'%(outputDirectory),'wb') )
    print('created dir: {o} '.format(o= outputDirectory))
    print('bining created successfully... ')
    print('Note than any additional call to createBins will overwrite directory {o}'.format(o= outputDirectory))
    sys.exit(0)

tnpBins = pickle.load( open( '%s/bining.pkl'%(outputDirectory),'rb') )


####################################################################
##### Create Histograms
####################################################################
for s in samplesDef.keys():
    sample =  samplesDef[s]
    if sample is None: continue
    setattr( sample, 'tree'     ,'%s/fitter_tree' % tnpTreeDir )
    setattr( sample, 'histFile' , '%s/%s_%s.root' % ( outputDirectory , sample.name, options.flag ) )


if options.createHists:
## parallel    for sampleType in samplesDef.keys():
## parallel        sample =  samplesDef[sampleType]
## parallel        if sample is None : continue
## parallel        if sampleType == options.sample or options.sample == 'all' :
## parallel            print('creating histogram for sample ')
## parallel            sample.dump()
## parallel            var = { 'name' : 'pair_mass', 'nbins' : 60, 'min' : 60, 'max': 120 }
## parallel            if sample.mcTruth:
## parallel                var = { 'name' : 'pair_mass', 'nbins' : 60, 'min' : 60, 'max': 120 }
## parallel            tnpRoot.makePassFailHistograms( sample, flags[options.flag][0], tnpBins['bins'], binningDef, flags[options.flag][1], var)#tnpBins, var )
    def parallel_hists(sampleType):
        # if not sampleType == 'mcNom': return
        sample = samplesDef[sampleType]
        if sample is not None and (sampleType == options.sample or options.sample == 'all'):
            print('creating histogram for sample ')
            sample.dump()
            var = { 'name' : 'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            if sample.mcTruth:
                var = { 'name' : 'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            ## tnpRoot.makePassFailHistograms( sample, tnpConf.flags[args.flag], tnpBins, var )
            tnpRoot.makePassFailHistograms( sample, flags[options.flag][0], tnpBins['bins'], binningDef, flags[options.flag][1], var)#tnpBins, var )
    
    pool = Pool()
    pool.map(parallel_hists, samplesDef.keys())

    sys.exit(0)


####################################################################
##### Actual Fitter
####################################################################
sampleToFit = samplesDef['data']
if sampleToFit is None:
    print('[tnpEGM_fitter, prelim checks]: sample (data or MC) not available... check your settings')
    sys.exit(1)

sampleMC = samplesDef['mcNom']

if sampleMC is None:
    print('[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings')
    sys.exit(1)
for s in samplesDef.keys():
    sample =  samplesDef[s]
    if sample is None: continue
    setattr( sample, 'mcRef'     , sampleMC )
    setattr( sample, 'nominalFit', '%s/%s_%s.nominalFit.root' % ( outputDirectory , sample.name, options.flag ) )
    setattr( sample, 'altSigFit' , '%s/%s_%s.altSigFit.root'  % ( outputDirectory , sample.name, options.flag ) )
    setattr( sample, 'altBkgFit' , '%s/%s_%s.altBkgFit.root'  % ( outputDirectory , sample.name, options.flag ) )



### change the sample to fit is mc fit
if options.mcSig :
    sampleToFit = samplesDef['mcNom']

if  options.doFit:
    print('sampleToFit.dump()', sampleToFit.dump())
    sampleToFit.dump()
    ## parallel for ib in range(len(tnpBins['bins'])):
    def parallel_fit(ib): ## parallel
        if (options.binNumber >= 0 and ib == options.binNumber) or options.binNumber < 0:
            if options.altSig:                 
                fitUtils.histFitterAltSig(  sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit )
            elif options.altBkg:
                fitUtils.histFitterAltBkg(  sampleToFit, tnpBins['bins'][ib], tnpParAltBkgFit )
            else:
                fitUtils.histFitterNominal( sampleToFit, tnpBins['bins'][ib], tnpParNomFit )

    pool = Pool() ## parallel
    pool.map(parallel_fit, range(len(tnpBins['bins']))) ## parallel


    options.doPlot = True
     
####################################################################
##### dumping plots
####################################################################
if  options.doPlot:
    fileName = sampleToFit.nominalFit
    fitType  = 'nominalFit'
    if options.altSig : 
        fileName = sampleToFit.altSigFit
        fitType  = 'altSigFit'
    if options.altBkg : 
        fileName = sampleToFit.altBkgFit
        fitType  = 'altBkgFit'

    os.system('hadd -f %s %s' % (fileName, fileName+'_bin_bin*')) #fileName.replace('.root', '*.root')))
    os.system('sleep 3')
    os.system('rm '+fileName+'_bin_bin*')
        
    plottingDir = '%s/plots/%s/%s' % (outputDirectory,sampleToFit.name,fitType)
    if not os.path.exists( plottingDir ):
        os.makedirs( plottingDir )
    shutil.copy('etc/inputs/index.php.listPlots','%s/index.php' % plottingDir)

    for ib in range(len(tnpBins['bins'])):
        if (options.binNumber >= 0 and ib == options.binNumber) or options.binNumber < 0:
            tnpRoot.histPlotter( fileName, tnpBins['bins'][ib], plottingDir )

    print(' ===> Plots saved in <=======')
#    print('localhost/%s/' % plottingDir)


####################################################################
##### dumping egamma txt file 
####################################################################
if options.sumUp:
    sampleToFit.dump()
    info = {
        'data'        : sampleToFit.histFile,
        'dataNominal' : sampleToFit.nominalFit,
        'dataAltSig'  : sampleToFit.altSigFit ,
        'dataAltBkg'  : sampleToFit.altBkgFit ,
        'mcNominal'   : sampleToFit.mcRef.histFile,
        'mcAlt'       : None,
        'tagSel'      : None
        }

    if not samplesDef['mcAlt' ] is None:
        info['mcAlt'    ] = samplesDef['mcAlt' ].histFile
    if not samplesDef['tagSel'] is None:
        info['tagSel'   ] = samplesDef['tagSel'].histFile

    print(info)

    effis = None
    effFileName ='%s/egammaEffi.txt' % outputDirectory 
    fOut = open( effFileName,'w')
    
    for ib in range(len(tnpBins['bins'])):
        effis = tnpRoot.getAllEffi( info, tnpBins['bins'][ib] )

        ### formatting assuming 2D bining -- to be fixed        
        v1Range = tnpBins['bins'][ib]['title'].split(';')[1].split('<')
        v2Range = tnpBins['bins'][ib]['title'].split(';')[2].split('<')
        if ib == 0 :
            astr = '### var1 : %s' % v1Range[1]
            print(astr)
            fOut.write( astr + '\n' )
            astr = '### var2 : %s' % v2Range[1]
            print(astr)
            fOut.write( astr + '\n' )
            
        astr =  '%+8.3f\t%+8.3f\t%+8.3f\t%+8.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f' % (
            float(v1Range[0]), float(v1Range[2]),
            float(v2Range[0]), float(v2Range[2]),
            effis['dataNominal'][0],effis['dataNominal'][1],
            effis['mcNominal'  ][0],effis['mcNominal'  ][1],
            effis['dataAltBkg' ][0],
            effis['dataAltSig' ][0],
            effis['mcAlt' ][0],
            effis['tagSel'][0],
            )
        print(astr)
        fOut.write( astr + '\n' )
    fOut.close()

    print('Effis saved in file : ',  effFileName)
    import libPython.EGammaID_scaleFactors as egm_sf
    egm_sf.doEGM_SFs(effFileName,sampleToFit.lumi)
