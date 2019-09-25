#!/bin/env python
## USAGE (for fits): python etc/scripts/submitBootstrapHists.py -n 500 -r 4 -s doFit --outdir condor_out

jobstring = '''#!/bin/sh
ulimit -c 0 -S
ulimit -c 0 -H
set -e
export SCRAM_ARCH=slc6_amd64_gcc630
cd CMSSWBASE
eval `scramv1 runtime -sh`
cd WORKDIR
TNPSTRING

'''

def makeCondorFile(jobdir, srcFiles, options, logdir, errdir, outdirCondor):
    dummy_exec = open(jobdir+'/dummy_exec.sh','w')
    dummy_exec.write('#!/bin/bash\n')
    dummy_exec.write('bash $*\n')
    dummy_exec.close()
     
    condor_file_name = jobdir+'/condor_submit.condor'
    condor_file = open(condor_file_name,'w')
    condor_file.write('''Universe = vanilla
Executable = {de}
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = {ld}/$(ProcId).log
Output     = {od}/$(ProcId).out
Error      = {ed}/$(ProcId).error
getenv      = True
next_job_start_delay = 1
environment = "LS_SUBCWD={here}"
+MaxRuntime = {rt}\n
'''.format(de=os.path.abspath(dummy_exec.name), ld=os.path.abspath(logdir), od=os.path.abspath(outdirCondor),ed=os.path.abspath(errdir),
           rt=int(options.runtime*3600), here=os.environ['PWD'] ) )
    if os.environ['USER'] in ['mdunser', 'psilva']:
        condor_file.write('+AccountingGroup = "group_u_CMST3.all"\n\n\n')
    if options.nThreads:
        condor_file.write('request_cpus = {nt} \n\n'.format(nt=options.nThreads))
    for sf in srcFiles:
        condor_file.write('arguments = {sf} \nqueue 1 \n\n'.format(sf=os.path.abspath(sf)))
    condor_file.close()
    return condor_file_name


import ROOT, random, array, os, sys

if __name__ == "__main__":
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option(        '--nBins'      , dest='nBins'         , type=int           , default=24     , help='make 1 job / TnP bin (data is extremely slow)')
    parser.add_option('-n'  , '--nreplicas'  , dest='nReplicas'     , type=int           , default=1      , help='make 1 data and MC replica for each job')
    parser.add_option('-t'  , '--threads'    , dest='nThreads'      , type=int           , default=None   , help='use nThreads in the fit (suggested 2 for single charge, 1 for combination)')
    parser.add_option('-r'  , '--runtime'    , default=8            , type=int                            , help='New runtime for condor resubmission in hours. default None: will take the original one.');
    parser.add_option('-s'  , '--step'       , dest='step'          , type="string", default='createHists', help='step for the tnp tool');
    parser.add_option(        '--outdir'     , dest='outdir'        , type="string", default=None,          help='outdirectory');
    (options, args) = parser.parse_args()

    toolStep = '--'+options.step

    absopath  = os.path.abspath(options.outdir)
    if not options.outdir:
        raise RuntimeError, 'ERROR: give at least an output directory. there will be a HUGE number of jobs!'
    else:
        if not os.path.isdir(absopath):
            print 'making a directory and running in it'
            os.system('mkdir -p {od}'.format(od=absopath))

    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir {od}'.format(od=jobdir))
    logdir = absopath+'/logs/'
    if not os.path.isdir(logdir):
        os.system('mkdir {od}'.format(od=logdir))
    errdir = absopath+'/errs/'
    if not os.path.isdir(errdir):
        os.system('mkdir {od}'.format(od=errdir))
    outdirCondor = absopath+'/outs/'
    if not os.path.isdir(outdirCondor):
        os.system('mkdir {od}'.format(od=outdirCondor))

    srcfiles = []
    ijob = 0
    for j in xrange(int(options.nReplicas)):
        for b in range(int(options.nBins)):
            ## make new file for evert parameter and point
            job_file_name = jobdir+'/job_{ij}.sh'.format(ij=ijob)
            log_file_name = logdir+'/job_{ij}.log'.format(ij=ijob)
            tmp_file = open(job_file_name, 'w')
            
            tmp_filecont = jobstring
            cmd = 'python scaleEGM_fitter.py etc/config/settings_elScale_allEras.py --flag ScaleFullID {toolstep} --iResample {res} --iBin {bin} --batch '.format(cmssw=os.environ['CMSSW_BASE'],toolstep=toolStep,res=j,bin=b)
            tmp_filecont = tmp_filecont.replace('TNPSTRING', cmd)
            tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
            tmp_filecont = tmp_filecont.replace('WORKDIR', os.environ['CMSSW_BASE']+'/src/egm_tnp_analysis/')
            tmp_file.write(tmp_filecont)
            tmp_file.close()
            srcfiles.append(job_file_name)
            ijob += 1
    cf = makeCondorFile(jobdir,srcfiles,options, logdir, errdir, outdirCondor)
    subcmd = 'condor_submit {rf} '.format(rf = cf)

    print subcmd

    sys.exit()

