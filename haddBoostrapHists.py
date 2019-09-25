#!/bin/env python

import pickle
import argparse
import os
import libPython.binUtils  as tnpBiner

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='tnp EGM fitter')
    parser.add_argument('--inputDir'       , dest = 'inputDir'   ,              default='./' , help='Directory with the histograms to be hadded')
    parser.add_argument('--nResamples'     , dest = 'resamples'  ,  type= int,  default=500  , help='Number of bootstrapping resamples')
    args = parser.parse_args()

    indir = args.inputDir
    tnpBins = pickle.load( open( '{indir}/bining.pkl'.format(indir=indir),'rb') )
    
    allfiles = list(filter(lambda x: x.endswith('.root'), os.listdir(indir)))
    datasets = []
    for f in allfiles: 
        basename = f.split('_stat')[0]
        if basename not in datasets and '.root' not in basename: datasets.append(basename)
    print "==> Datasets are: ",datasets
    bins = [tnpBins['bins'][i]['name'] for i in xrange(len(tnpBins['bins']))]
    print "==> TnP bins = ", bins

    pwd = os.environ['PWD']
    for istat in range(args.resamples):
        print "MERGING boostrap resample # ",istat
        haddtxt = 'haddfile.sh'
        haddfile = open(haddtxt,'w')
        for ds in datasets:
            haddcmd = 'hadd -f {indir}/{dataset}_stat{ist}.root '.format(indir=indir,dataset=ds,ist=istat)
            haddcmd += ' '.join(['{indir}/{dataset}_stat{ist}_{binname}.root'.format(indir=indir,dataset=ds,ist=istat,binname=name) for name in bins])
            haddcmd += '\n'
            haddfile.write(haddcmd)
        haddfile.close()
        os.system('source {pwd}/{cmdfile}'.format(pwd=pwd,cmdfile=haddtxt))

