import os
from submitBootstrapHists import jobstring,makeCondorFile

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option(        '--outdir'     , dest='outdir'        , type="string", default=None,          help='outdirectory');
    (options, args) = parser.parse_args()
    options.runtime = 1
    options.nThreads = 1

    absopath  = os.path.abspath(options.outdir)
    if not options.outdir:
        raise RuntimeError, 'ERROR: give at least an output directory. there will be a HUGE number of jobs!'
    else:
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
    for j in range(64):
        job_file_name = jobdir+'/job_{ij}.sh'.format(ij=j)
        log_file_name = logdir+'/job_{ij}.log'.format(ij=j)
        tmp_file = open(job_file_name, 'w')
        
        tmp_filecont = jobstring
        cmd = 'python scaleEGM_fitter.py settings_tmp.py --flag ScaleHighPt --doFit --altSig --iBin {ibin} --batch '.format(ibin=j)
        tmp_filecont = tmp_filecont.replace('TNPSTRING', cmd)
        tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
        tmp_filecont = tmp_filecont.replace('WORKDIR', os.environ['CMSSW_BASE']+'/src/egm_tnp_analysis/')
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        srcfiles.append(job_file_name)

    cf = makeCondorFile(jobdir,srcfiles,options, logdir, errdir, outdirCondor)
    subcmd = 'condor_submit {rf} '.format(rf = cf)

    print(subcmd)
