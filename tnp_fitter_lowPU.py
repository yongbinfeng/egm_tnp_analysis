
### python specific import
import os
import sys
import pickle
import shutil
from multiprocessing import Pool
import datetime

def joinCuts(_list):
    return '(' + ' && '.join(_list) + ')'

from optparse import OptionParser
#parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
#parser.add_option('-c', '--channel'   , dest='channel'   , type='string'      , default='mu',  help='run tnp ntuples for muons/electrons. default mu')
#parser.add_option('-i', '--indir'     , dest='indir'     , type='string'      , default=''  ,  help='directory with the trees')


parser = OptionParser(usage="%prog [options] ")
parser.add_option('--checkBins'  , action='store_true'  , help = 'check  bining definition')
parser.add_option('--mcTruth'    , dest='mcTruth', default=False, action='store_true'  , help = 'run with mc truth matching on the leptons')
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

massbins, massmin, massmax = 80, 50, 130

## define the binning here, much easier...
binning_eta = [-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0.,0.3,0.9,1.2,1.6,2.1,2.4]
binning_pt  = [25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 65]#, 80]#, 10000]

if options.flag.startswith('el'):
    binning_eta = [-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4]
    binning_pt  = [25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 65]#, 80]#, 10000]
if not 'trigger' in options.flag:
    binning_pt  = [25, 35., 50., 65.]

## binning from yongbin
## Muon Channel:
## 12 eta bins: [-2.4,-2.1,-1.6,-1.2,-0.9,-0.3,0,0.3,0.9,1.2,1.6,2.1,2.4] 
## 12 pt bins for HLT: [25, 26.5, 28, 29.5, 31, 32.5, 35, 40, 45, 50, 60, 80, 10000]
## 3 pt bins for selections, etc: [25, 35, 50, 10000]
## 
## Electron Channel: 
## 12 eta bins: [-2.4,-2.0,-1.566,-1.4442,-1.0,-0.5,0,0.5,1.0,1.4442,1.566,2.0,2.4]
## Same pt bins as the muons for HLT and selections.


binningDef = {
    'eta' : {'var' : 'probe_eta', 'type': 'float', 'bins': binning_eta},
    'pt'  : {'var' : 'probe_pt' , 'type': 'float', 'bins': binning_pt }
}

## now the flags and the cuts, pretty straight forward
cutBase     = 'tag_pt > 25. && probe_pt > 25. && probe_isMuon > 0.5'

if options.mcTruth:
    cutBase += ' && tag_matchMC == 1 && probe_matchMC == 1'

cutMinus       = 'probe_charge < 0'
cutPlus        = 'probe_charge > 0'


cutTrigger     = 'probe_triggerMatch > 0.5'
cutIdIp        = 'probe_mediumId > 0.5 && fabs(probe_dxy) < 0.05 && probe_isGlobal > 0.5'
cutIso         = 'probe_iso < 0.15'

cutIdIso       = 'probe_mediumId >= 3'
## any additional cut to apply
additionalCuts = None

## the first item in the flags is the flag to test, the second is the base cut to be applied to all events
flags = {
    'mu_idip_both'      : (cutIdIp     , joinCuts([cutBase          ]) ),
    'mu_idip_minus'     : (cutIdIp     , joinCuts([cutBase, cutMinus]) ),
    'mu_idip_plus'      : (cutIdIp     , joinCuts([cutBase, cutPlus ]) ),

    'mu_iso_both'       : (cutIso      , joinCuts([cutBase          , cutIdIp]) ),
    'mu_iso_minus'      : (cutIso      , joinCuts([cutBase, cutMinus, cutIdIp]) ),
    'mu_iso_plus'       : (cutIso      , joinCuts([cutBase, cutPlus , cutIdIp]) ),

    'mu_trigger_both'   : (cutTrigger  , joinCuts([cutBase          , cutIdIp, cutIso]) ),
    'mu_trigger_minus'  : (cutTrigger  , joinCuts([cutBase, cutMinus, cutIdIp, cutIso]) ),
    'mu_trigger_plus'   : (cutTrigger  , joinCuts([cutBase, cutPlus , cutIdIp, cutIso]) ),

    'el_idiso_both'     : (cutIdIso     , joinCuts([cutBase          ]) ),
    'el_idiso_minus'    : (cutIdIso     , joinCuts([cutBase, cutMinus]) ),
    'el_idiso_plus'     : (cutIdIso     , joinCuts([cutBase, cutPlus ]) ),

    'el_trigger_both'   : (cutTrigger  , joinCuts([cutBase          , cutIdIso]) ),
    'el_trigger_minus'  : (cutTrigger  , joinCuts([cutBase, cutMinus, cutIdIso]) ),
    'el_trigger_plus'   : (cutTrigger  , joinCuts([cutBase, cutPlus , cutIdIso]) )

}


from libPython.tnpClassUtils import tnpSample
inputDir = '/data/shared/tnp/lowPU/2021-10-19-firstGo/'
wmass_selection = {
    'mu_DY_lowPU'    : tnpSample('mu_DY_lowPU'    , inputDir +'/dy_mumu.root'    , isMC = True       , nEvts = -1 ) ,
    'mu_data_lowPU'  : tnpSample('mu_data_lowPU'  , inputDir +'/data_mumu.root'  , lumi =  0.200 ) ,
    'el_DY_lowPU'    : tnpSample('el_DY_lowPU'    , inputDir +'/dy_elel.root'    , isMC = True       , nEvts = -1 ) ,
    'el_data_lowPU'  : tnpSample('el_data_lowPU'  , inputDir +'/data_elel.root'  , lumi =  0.200 ) ,

}

samples_data_mumu  = wmass_selection['mu_data_lowPU'  ].clone()
samples_dy_mumu    = wmass_selection['mu_DY_lowPU' ].clone()
samples_data_elel  = wmass_selection['el_data_lowPU'  ].clone()
samples_dy_elel    = wmass_selection['el_DY_lowPU' ].clone()

samples_dy   = samples_dy_mumu
samples_data = samples_data_mumu
if options.flag.startswith('el'):
    samples_dy   = samples_dy_elel
    samples_data = samples_data_elel

weightName   = 'std::copysign(1.0f,genWeight)'

#weightName = 'totWeight/std::abs(totWeight)*std::max(4.f,std::abs(totWeight))'


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
#    "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
#    "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
#    "acmsP[60.,50.,80.]","betaP[0.05,0.01,0.08]","gammaP[0.1, 0, 1]","peakP[90.0]",
#    "acmsF[60.,50.,80.]","betaF[0.05,0.01,0.08]","gammaF[0.1, 0, 1]","peakF[90.0]",
    "meanP[-0.0,-5.0,5.0]","sigmaP[0.5,0.1,5.0]",
    "meanF[-0.0,-5.0,5.0]","sigmaF[0.5,0.1,5.0]",
    "acmsP[60.,40.,130.]","betaP[0.05,0.01,0.08]","gammaP[0.1, 0, 1]","peakP[90.0]",
    "acmsF[60.,40.,130.]","betaF[0.05,0.01,0.08]","gammaF[0.1, 0, 1]","peakF[90.0]",
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

today = datetime.date.isoformat(datetime.date.today())

baseOutDir = 'resultsLowPU_nodz_dxybs{mc}_{d}/efficiencies_{era}/'.format(d=today,era=options.era,mc='_mcTruth' if options.mcTruth else '')

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
    setattr( sample, 'tree'     ,'fitter_tree' )
    setattr( sample, 'histFile' , '%s/%s_%s.root' % ( outputDirectory , sample.name, options.flag ) )


if options.createHists:
    def parallel_hists(sampleType):
        # if not sampleType == 'mcNom': return
        sample = samplesDef[sampleType]
        if sample is not None and (sampleType == options.sample or options.sample == 'all'):
            print('creating histogram for sample ')
            sample.dump()
            if 'tracking' in options.flag:
                var = { 'namePassing' : 'pair_mass','nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            elif 'alttrack' in options.flag:
                var = { 'namePassing' : 'pair_massMatchedTrack','nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            else:
                var = { 'name' : 'pair_mass', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
            if sample.mcTruth:
                if 'tracking' in options.flag:
                    var = { 'namePassing' : 'pair_mass','nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
                elif 'alttrack' in options.flag:
                    var = { 'namePassing' : 'pair_massMatchedTrack','nameFailing' : 'pair_massStandalone', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
                else:
                    var = { 'name' : 'pair_mass', 'nbins' : massbins, 'min' : massmin, 'max': massmax }
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
    usePassingTemplateForNominalFit = False
    sampleToFit.dump()
    ## parallel for ib in range(len(tnpBins['bins'])):
    def parallel_fit(ib): ## parallel
        if (options.binNumber >= 0 and ib == options.binNumber) or options.binNumber < 0:
            if options.altSig:                 
                fitUtils.histFitterAltSig(  sampleToFit, tnpBins['bins'][ib], tnpParAltSigFit, massbins, massmin, massmax )
            elif options.altBkg:
                fitUtils.histFitterAltBkg(  sampleToFit, tnpBins['bins'][ib], tnpParAltBkgFit, massbins, massmin, massmax )
            else:
                fitUtils.histFitterNominal( sampleToFit, tnpBins['bins'][ib], tnpParNomFit, massbins, massmin, massmax, usePassingTemplateForNominalFit )

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
    os.system('sleep 1')
    os.system('rm '+fileName+'_bin_bin*')
        
    plottingDir = '%s/plots/%s/%s' % (outputDirectory,sampleToFit.name,fitType)
    if not os.path.exists( plottingDir ):
        os.makedirs( plottingDir )
    shutil.copy('etc/inputs/index.php.listPlots','%s/index.php' % plottingDir)

    verbosePlotting = False
    for ib in range(len(tnpBins['bins'])):
        if (options.binNumber >= 0 and ib == options.binNumber) or options.binNumber < 0:
            tnpRoot.histPlotter( fileName, tnpBins['bins'][ib], plottingDir, -1, verbosePlotting ) ## the -1 is form marc, something with replica

    print(' ===> Plots saved in <=======')
#    print('localhost/%s/' % plottingDir)


####################################################################
##### dumping egamma txt file 
####################################################################
if options.sumUp:
    import ROOT, copy
    from pprint import pprint
    #print('this is the dump of sampleToFit:')
    sampleToFit.dump()
    #pprint(vars(sampleToFit.mcRef))
    #print('done with dump')
    info = {
        'data'        : sampleToFit.histFile,
        'dataNominal' : sampleToFit.nominalFit,
        'dataAltSig'  : sampleToFit.altSigFit ,
        'dataAltBkg'  : sampleToFit.altBkgFit ,
        'mcNominal'   : sampleToFit.mcRef.histFile,
        ## marc 'mcAlt'       : None,
        'mcAlt'       : sampleToFit.mcRef.altSigFit,
        'tagSel'      : None
        }

    if not samplesDef['mcAlt' ] is None:
        info['mcAlt'    ] = samplesDef['mcAlt' ].histFile
    if not samplesDef['tagSel'] is None:
        info['tagSel'   ] = samplesDef['tagSel'].histFile

    print(info)

    effis = None
    effFileName ='%s/allEfficiencies.txt' % outputDirectory 
    fOut = open( effFileName,'w')
    
    for ib,_bin in enumerate(tnpBins['bins']):
        effis = tnpRoot.getAllEffi( info, _bin )

        ### formatting assuming 2D bining -- to be fixed        
        v1Range = _bin['title'].split(';')[1].split('<')
        v2Range = _bin['title'].split(';')[2].split('<')
        
        if not ib:
            astr = '### var1 : %s\n' % v1Range[1]
            fOut.write( astr )
            astr = '### var2 : %s\n' % v2Range[1]
            fOut.write( astr )
            exp = '{v0:8s}\t{v1:8s}\t{v2:8s}\t{v3:8s}\t{edv:10s}\t{ede:10s}\t{emcv:10s}\t{emce:10s}\t{edalts:15s}\t{edaltse:15s}\t{emcalt:15s}\t{emcalte:15s}\t{edaltb:15s}\t{etagsel:10s}\n'.format(
                v0='var1min', v1='var1max', v2='var2min', v3='var2max', 
                edv='eff data', ede='err data', 
                emcv='eff mc', emce='err mc', 
                edalts='eff data altS', edaltse='err data altS',
                emcalt='eff mc alt', emcalte='err mc alt', edaltb='eff data altB', etagsel='eff tag sel')
            print(exp)
            fOut.write(exp)
            
        astr =  '%-+8.3f\t%-+8.3f\t%-+8.3f\t%-+8.3f\t%-10.5f\t%-10.5f\t%-10.5f\t%-10.5f\t%-15.5f\t%-15.5f\t%-15.5f\t%-15.5f\t%-15.5f\t%-10.5f' % (
            float(v1Range[0]), float(v1Range[2]),
            float(v2Range[0]), float(v2Range[2]),
            effis['dataNominal'][0],effis['dataNominal'][1],
            effis['mcNominal'  ][0],effis['mcNominal'  ][1],
            effis['dataAltSig' ][0],effis['dataAltSig' ][1],
            effis['mcAlt' ][0], effis['mcAlt' ][1],
            effis['dataAltBkg' ][0],
            effis['tagSel'][0],
            )
        print(astr)
        fOut.write( astr + '\n' )

        canv_all = ROOT.TCanvas(_bin['name'], _bin['name'], 1200, 1200)
        canv_all.Divide(3,3)
        canv_all.Draw()
        ipad = 1
        canv_all.cd(0)
        txt = ROOT.TLatex(); txt.SetTextFont(42); txt.SetTextSize(0.03)
        txt.SetNDC()
        txt.DrawLatex(0.01, 0.97, '{n}'.format(n=_bin['name'].replace('_',' ').replace('To', '-').replace('probe ', '').replace('m','-').replace('pt','XX').replace('p','.').replace('XX','p_{T}')))
        txt.SetTextSize(0.07)
        for ip, p in enumerate(effis['canv_mcAlt'].GetListOfPrimitives()):
            if not ip: continue
            canv_all.cd(ipad)
            p.SetPad(0.05, 0.00, 0.95, 0.90)
            p.Draw()
            ipad+=1
        canv_all.cd(ipad)
        txt.DrawLatex(0.01, 0.85, 'MC counting efficiency:')
        tmp = effis['mcNominal']
        txt.DrawLatex(0.20, 0.75, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
        txt.DrawLatex(0.20, 0.65, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.55, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
        txt.SetTextFont(42)
        tmp = effis['mcAlt']
        txt.DrawLatex(0.01, 0.35, 'MC fitted signal:')
        txt.DrawLatex(0.20, 0.25, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
        txt.DrawLatex(0.20, 0.15, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.05, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
        txt.SetTextFont(42)
        ipad+=1
        for ip, p in enumerate(effis['canv_dataNominal'].GetListOfPrimitives()):
            if not ip: continue
            canv_all.cd(ipad)
            p.SetPad(0.05, 0.00, 0.95, 0.90)
            p.Draw()
            ipad+=1
        canv_all.cd(ipad)
        tmp = effis['dataNominal']
        txt.DrawLatex(0.01, 0.65, 'data nominal:')
        txt.DrawLatex(0.20, 0.55, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
        txt.DrawLatex(0.20, 0.45, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.35, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
        txt.SetTextFont(42)
        ipad+=1
        for ip, p in enumerate(effis['canv_dataAltSig'].GetListOfPrimitives()):
            if not ip: continue
            canv_all.cd(ipad)
            p.SetPad(0.05, 0.00, 0.95, 0.90)
            p.Draw()
            ipad+=1
        canv_all.cd(ipad)
        tmp = effis['dataAltSig']
        txt.DrawLatex(0.01, 0.65, 'data alternative:')
        txt.DrawLatex(0.20, 0.55, 'passing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[2],ne=tmp[4]))
        txt.DrawLatex(0.20, 0.45, 'failing: {n:.1f} #pm {ne:.1f}'.format(n=tmp[3],ne=tmp[5]))
        txt.SetTextFont(62)
        txt.DrawLatex(0.20, 0.35, 'efficiency: {e:.2f} #pm {ee:.2f} %'.format(e=tmp[0]*100., ee=tmp[1]*100.))
        txt.SetTextFont(42)

        #effis['canv_dataAltSig'].Draw()

        odllevel = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        canv_all.SaveAs(outputDirectory+'/plots/{n}_all.pdf'.format(n=_bin['name']))
        canv_all.SaveAs(outputDirectory+'/plots/{n}_all.png'.format(n=_bin['name']))
        ROOT.gErrorIgnoreLevel = odllevel
    os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {d}/index.php'.format(d=outputDirectory+'/plots/'))
        

    fOut.close()

    print('Efficiencies saved in file : ',  effFileName)
    import libPython.EGammaID_scaleFactors as makesf
    makesf.doSFs(effFileName,sampleToFit.lumi,['pt', 'eta'], outputDirectory+'/plots/')
    ## put here the new file for marco
