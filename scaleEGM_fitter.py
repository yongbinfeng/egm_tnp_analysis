## steps similar to TnP, but the histogram making is multiplied by nResamples (i.e. 500), so this step has to be submitted in batch.
# python scaleEGM_fitter.py etc/config/settings_elScale_allEras.py --flag ScaleFullID --checkBins
# python scaleEGM_fitter.py etc/config/settings_elScale_allEras.py --flag ScaleFullID --createBins
# python etc/scripts/submitBootstrapHists.py -n 500 -r 8 --outdir condor_out # (this runs the --createHists step)
# rest is still not implemented......... 

### python specific import
import argparse
import os
import sys
import pickle
import shutil
import json

parser = argparse.ArgumentParser(description='tnp EGM fitter')
parser.add_argument('--checkBins'  , action='store_true'  , help = 'check  bining definition')
parser.add_argument('--createBins' , action='store_true'  , help = 'create bining definition')
parser.add_argument('--createHists', action='store_true'  , help = 'create histograms')
parser.add_argument('--sample'     , default='all'        , help = 'create histograms (per sample, expert only)')
parser.add_argument('--altSig'     , action='store_true'  , help = 'alternate signal model fit')
parser.add_argument('--altBkg'     , action='store_true'  , help = 'alternate background model fit')
parser.add_argument('--doFit'      , action='store_true'  , help = 'fit sample (sample should be defined in settings.py)')
parser.add_argument('--mcSig'      , action='store_true'  , help = 'fit MC nom [to init fit parama]')
parser.add_argument('--doPlot'     , action='store_true'  , help = 'plotting')
parser.add_argument('--sumUp'      , action='store_true'  , help = 'sum up efficiencies')
parser.add_argument('--iBin'       , dest = 'binNumber'   , type = int,  default=-1, help='bin number (to refit individual bin)')
parser.add_argument('--flag'       , default = None       , help ='WP to test')
parser.add_argument('settings'     , default = None       , help = 'setting file [mandatory]')
parser.add_argument('--iResample'  , dest = 'statResample' , type = int,  default=-1, help='resample number for bootstrapping')
parser.add_argument('--batch'      , action='store_true'  , help = 'send in batch (it makes one root file/bin)')


args = parser.parse_args()

print '===> settings %s <===' % args.settings
importSetting = 'import %s as tnpConf' % args.settings.replace('/','.').split('.py')[0]
print importSetting
exec(importSetting)

### tnp library
import libPython.binUtils  as tnpBiner
import libPython.rootUtils as tnpRoot


if args.flag is None:
    print '[tnpEGM_fitter] flag is MANDATORY, this is the working point as defined in the settings.py'
    sys.exit(0)
    
if not args.flag in tnpConf.flags.keys() :
    print '[tnpEGM_fitter] flag %s not found in flags definitions' % args.flag
    print '  --> define in settings first'
    print '  In settings I found flags: '
    print tnpConf.flags.keys()
    sys.exit(1)

outputDirectory = '%s/%s/' % (tnpConf.baseOutDir,args.flag)

print '===>  Output directory: '
print outputDirectory


####################################################################
##### Create (check) Bins
####################################################################
if args.checkBins:
    tnpBins = tnpBiner.createBins(tnpConf.biningDef,tnpConf.cutBase)
    tnpBiner.tuneCuts( tnpBins, tnpConf.additionalCuts )
    for ib in range(len(tnpBins['bins'])):
        print tnpBins['bins'][ib]['name']
        print '  - cut: ',tnpBins['bins'][ib]['cut']
    sys.exit(0)
    
if args.createBins:
    if os.path.exists( outputDirectory ):
            shutil.rmtree( outputDirectory )
    os.makedirs( outputDirectory )
    tnpBins = tnpBiner.createBins(tnpConf.biningDef,tnpConf.cutBase)
    tnpBiner.tuneCuts( tnpBins, tnpConf.additionalCuts )
    pickle.dump( tnpBins, open( '%s/bining.pkl'%(outputDirectory),'wb') )
    print 'created dir: %s ' % outputDirectory
    print 'bining created successfully... '
    print 'Note than any additional call to createBins will overwrite directory %s' % outputDirectory
    sys.exit(0)

tnpBins = pickle.load( open( '%s/bining.pkl'%(outputDirectory),'rb') )


####################################################################
##### Create Histograms
####################################################################
# put here all the replicas which have histograms for MC and data (because of failed jobs, this could be < tnpConf.nResamples
goodReplicas = []
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    print '====> mkdir -p {odir}/{sample}'.format(odir=outputDirectory, sample=sample.name)
    os.system('mkdir -p {odir}/{sample}'.format(odir=outputDirectory, sample=sample.name))
    setattr( sample, 'tree'     ,'%s/fitter_tree' % tnpConf.tnpTreeDir )
    for ir in xrange(tnpConf.nResamples):
        fname = '%s/%s/%s_%s_stat%d.root' % ( outputDirectory, sample.name, sample.name, args.flag, ir )
        checkHist = True
        if not args.createHists:
            if not os.path.exists(fname):
                print "WARNING! ",fname, " does not exist. Skipping this resample."
                checkHist = False
        if checkHist:
            setattr( sample, 'histFile%d' % ir , fname )

# put here all the replicas which have histograms for MC and data (because of failed jobs, this could be < tnpConf.nResamples
goodReplicas = range(tnpConf.nResamples)
fileWithGoodReplicas = 'goodReplicas.txt'
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    for ir in xrange(tnpConf.nResamples):
        if ir in goodReplicas and not hasattr(sample, 'histFile%d' % ir): goodReplicas.remove(ir)
print "There are ", len(goodReplicas), " good replicas = ", goodReplicas
print "Writing good replicas into ",fileWithGoodReplicas
fileReplicas = open(fileWithGoodReplicas,'w')
fileReplicas.write(json.dumps(goodReplicas))
fileReplicas.close()

if args.createHists:
    if (args.binNumber >= 0 and args.batch):
        tnpBins['bins'] = [tnpBins['bins'][args.binNumber]]
        print "Making histograms for the unique bin: ",tnpBins['bins']
    for sampleType in tnpConf.samplesDef.keys():
        sample =  tnpConf.samplesDef[sampleType]
        if sample is None : continue
        if sampleType == args.sample or args.sample == 'all' :
            print 'creating histogram for sample '
            sample.dump()
            var = { 'name' : 'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            if sample.mcTruth:
                var = { 'name' : 'pair_mass', 'nbins' : 80, 'min' : 50, 'max': 130 }
            if args.batch and args.binNumber >= 0: 
                fname = getattr(sample,'histFile%d' % args.statResample)
                fname = fname.replace('.root','_%s.root' % tnpBins['bins'][0]['name'])
                setattr(sample,'histFile%d' % args.statResample, fname)
            tnpRoot.makeBootstrapHistograms( sample, tnpConf.flags[args.flag], tnpBins, var, args.statResample)

    sys.exit(0)


####################################################################
##### Actual Fitter
####################################################################
sampleToFit = tnpConf.samplesDef['data']
if sampleToFit is None:
    print '[tnpEGM_fitter, prelim checks]: sample (data or MC) not available... check your settings'
    sys.exit(1)

sampleMC = tnpConf.samplesDef['mcNom']

refReplica = goodReplicas[0]

if sampleMC is None:
    print '[tnpEGM_fitter, prelim checks]: MC sample not available... check your settings'
    sys.exit(1)
for s in tnpConf.samplesDef.keys():
    sample =  tnpConf.samplesDef[s]
    if sample is None: continue
    setattr( sample, 'mcRef'     , sampleMC )
    setattr( sample, 'nominalFit', '%s/%s_%s.nominalFit.root' % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altSigFit' , '%s/%s_%s.altSigFit.root'  % ( outputDirectory , sample.name, args.flag ) )
    setattr( sample, 'altBkgFit' , '%s/%s_%s.altBkgFit.root'  % ( outputDirectory , sample.name, args.flag ) )


### change the sample to fit is mc fit
if args.mcSig :
    sampleToFit = tnpConf.samplesDef['mcNom']

if  args.doFit:
    sampleToFit.dump()
    for ib in range(len(tnpBins['bins'])):
        if (args.binNumber >= 0 and ib == args.binNumber) or args.binNumber < 0:
            if args.altSig:                 
                tnpRoot.histScaleFitterAltSig(  sampleToFit, tnpBins['bins'][ib], tnpConf.tnpParAltSigFit, goodReplicas[0], batch=True ) # this fit takes long, may split by iBin
            elif args.altBkg:
                tnpRoot.histScaleFitterAltBkg(  sampleToFit, tnpBins['bins'][ib], tnpConf.tnpParAltBkgFit, goodReplicas[0] )
            else:
                if args.statResample>=0:
                    goodReplicas = [args.statResample] if args.statResample in goodReplicas else []
                for ir in goodReplicas:
                    if hasattr(sampleToFit,'histFile%d' % ir) and hasattr(sampleToFit.mcRef,'histFile%d' % ir):
                        print "FITTING replica ",ir
                        tnpRoot.histScaleFitterNominal( sampleToFit, tnpBins['bins'][ib], tnpConf.tnpParNomFit, ir, batch=True )
                    else: 
                        print "Replica ", ir, " skipped because it is missing either the data or MC replica"

    args.doPlot = True
     


####################################################################
##### dumping plots
####################################################################
if  args.doPlot:
    fileName = sampleToFit.nominalFit
    fitType  = 'nominalFit'
    if args.altSig : 
        fileName = sampleToFit.altSigFit
        fitType  = 'altSigFit'
    if args.altBkg : 
        fileName = sampleToFit.altBkgFit
        fitType  = 'altBkgFit'

    plottingDir = '%s/plots/%s/%s' % (outputDirectory,sampleToFit.name,fitType)
    if not os.path.exists( plottingDir ):
        os.makedirs( plottingDir )
    shutil.copy('etc/inputs/index.php.listPlots','%s/index.php' % plottingDir)

    for ib in range(len(tnpBins['bins'])):
        if (args.binNumber >= 0 and ib == args.binNumber) or args.binNumber < 0:
            for ir in goodReplicas:
                print "plotting replica = ", ir
                irFileName = fileName.replace('.root', '_Stat%d.root' % ir)
                tnpRoot.histPlotter( irFileName, tnpBins['bins'][ib], plottingDir, ir )

    print ' ===> Plots saved in <======='
    print '%s/' % plottingDir


####################################################################
##### dumping egamma txt file 
####################################################################
if args.sumUp:
    setattr( sampleToFit, 'nominalFit', '%s/%s_%s.nominalFit_Stat%d.root' % ( outputDirectory , sampleToFit.name, args.flag, refReplica ) )
    sampleToFit.dump()
    info = {
        'dataNominal' : sampleToFit.nominalFit,
        'dataAltSig'  : sampleToFit.altSigFit ,
        'dataAltBkg'  : sampleToFit.altBkgFit ,
        'mcAlt'       : None,
        'tagSel'      : None
        }

    for ir in goodReplicas:
        info['dataReplica%d' % ir] = '%s/%s_%s.nominalFit_Stat%d.root' % ( outputDirectory , sampleToFit.name, args.flag, ir )

    if not tnpConf.samplesDef['mcAlt' ] is None:
        info['mcAlt'    ] = tnpConf.samplesDef['mcAlt' ].histFile
    if not tnpConf.samplesDef['tagSel'] is None:
        info['tagSel'   ] = tnpConf.samplesDef['tagSel'].histFile

    scales = None
    scaleFileName ='%s/egammaScale.txt' % outputDirectory 
    fOut = open( scaleFileName,'w')
    
    for ib in range(len(tnpBins['bins'])):
        scales = tnpRoot.getAllScales( info, tnpBins['bins'][ib], refReplica )
        
        ### formatting assuming 2D bining -- to be fixed        
        v1Range = tnpBins['bins'][ib]['title'].split(';')[1].split('<')
        v2Range = tnpBins['bins'][ib]['title'].split(';')[2].split('<')
        if ib == 0 :
            astr = '### var1 : %s' % v1Range[1]
            print astr
            fOut.write( astr + '\n' )
            astr = '### var2 : %s' % v2Range[1]
            print astr
            fOut.write( astr + '\n' )
            
        astr =  '%+8.3f\t%+8.3f\t%+8.3f\t%+8.3f\t%5.3f\t%5.3f\t%5.3f' % (
            float(v1Range[0]), float(v1Range[2]),
            float(v2Range[0]), float(v2Range[2]),
            scales['dataNominal'][0],
            -1*scales['dataAltSig'][0], # the difinition of the bias is opposite
            scales['dataAltBkg'][0],
            )
        astr += '\t'+'\t'.join(['%+5.3f' % scales['dataReplica%d' % ir][0] for ir in goodReplicas])
        print astr
        fOut.write( astr + '\n' )
    fOut.close()

    print 'Effis saved in file : ',  scaleFileName
    import libPython.EGammaID_scaleSystematics as egm_scales
    egm_scales.doEGM_Scales(scaleFileName,sampleToFit.lumi)
