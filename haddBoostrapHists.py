#!/bin/env python

import pickle
import argparse
import os
import libPython.binUtils  as tnpBiner
from etc.scripts.submitBootstrapHists import jobstring,makeCondorFile

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='tnp EGM fitter')
    parser.add_argument('--runCheck'       , dest = 'runCheck'   ,  action='store_true',            default=False , help='Do not hadd, only check the existence of all the files')
    parser.add_argument('--inputDir'       , dest = 'inputDir'   ,              default='./' , help='Directory with the histograms to be hadded')
    parser.add_argument('--nResamples'     , dest = 'resamples'  ,  type= int,  default=100  , help='Number of bootstrapping resamples')
    parser.add_argument('--runtime'        , default=4           , type=int                            , help='New runtime for condor resubmission in hours. default None: will take the original one.');
    parser.add_argument('--threads'        , dest='nThreads'     , type=int           , default=None   , help='use nThreads in the fit (suggested 2 for single charge, 1 for combination)')
    parser.add_argument('--flag'           , dest='flag'         ,                    default='ScaleHighPt', help='Flag for the Scale');
    args = parser.parse_args()


    indir = args.inputDir
    tnpBins = pickle.load( open( '{indir}/bining.pkl'.format(indir=indir),'rb') )
    
    datasets = list(filter(lambda x: x!='plots' and not x.endswith('.pkl') and not x.endswith('.root'), os.listdir(indir)))
    print("==> Datasets are: ",datasets)
    bins = [tnpBins['bins'][i]['name'] for i in xrange(len(tnpBins['bins']))]
    #print("==> TnP bins = ", bins)

    if args.runCheck:
        cmds = []
        for istat in range(args.resamples):
            print("CHECKING boostrap resample # ",istat)
            for ds in datasets:
                for name in bins:
                    ib = name.split('_')[0].split('bin')[-1]
                    filename = '{indir}/{dataset}/{dataset}_{flag}_stat{ist}_{binname}.root'.format(indir=indir,dataset=ds,flag=args.flag,ist=istat,binname=name)
                    if not os.path.exists(filename):
                        cmd = 'python scaleEGM_fitter.py etc/config/settings_elScale_allEras.py --flag {flag} --createHists --iResample {ist} --iBin {ib} --batch'.format(ist=istat,flag=args.flag,ib=ib)
                        if cmd not in cmds: cmds.append(cmd)
        print("There are ",len(cmds)," hists to be redone.")

        absopath  = os.path.abspath('condor_hist_recover')        
        if not os.path.isdir(absopath):
            print('making a directory and running in it')
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
        for j in range(len(cmds)):
            job_file_name = jobdir+'/job_{ij}.sh'.format(ij=j)
            log_file_name = logdir+'/job_{ij}.log'.format(ij=j)
            tmp_file = open(job_file_name, 'w')
            
            tmp_filecont = jobstring
            tmp_filecont = tmp_filecont.replace('TNPSTRING', cmds[j])
            tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
            tmp_filecont = tmp_filecont.replace('WORKDIR', os.environ['CMSSW_BASE']+'/src/egm_tnp_analysis/')
            tmp_file.write(tmp_filecont)
            tmp_file.close()
            srcfiles.append(job_file_name)

        cf = makeCondorFile(jobdir,srcfiles,args, logdir, errdir, outdirCondor)
        subcmd = 'condor_submit {rf} '.format(rf = cf)
        print(subcmd)
        if len(srcfiles)==0:
            os.system('rm -r {od}'.format(od=jobdir))
            print("ALL FILES ARE GOOD.")
        exit(0)
        
    pwd = os.environ['PWD']
    for istat in range(args.resamples):
        print("MERGING boostrap resample # ",istat)
        haddtxt = 'haddfile.sh'
        haddfile = open(haddtxt,'w')
        for ds in datasets:
            haddcmd = 'hadd -f {indir}/{dataset}_{flag}_stat{ist}.root '.format(indir=indir,dataset=ds,flag=args.flag,ist=istat)
            haddcmd += ' '.join(['{indir}/{dataset}/{dataset}_{flag}_stat{ist}_{binname}.root'.format(indir=indir,dataset=ds,flag=args.flag,ist=istat,binname=name) for name in bins])
            haddcmd += '\n'
            haddfile.write(haddcmd)
        haddfile.close()
        os.system('source {pwd}/{cmdfile}'.format(pwd=pwd,cmdfile=haddtxt))

