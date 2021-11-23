import subprocess

def runCommands(flav,t,e,c,mc=False,lowPU=False):
    opt_e = '--era='+e
    opt_f = '--flag={f}_{t}_{c}'.format(f=flav,t=t,c=c)
    opt_m = '--mcTruth' if mc else ''
    cmds = []
    ex = 'tnpEGM_fitter.py' if not lowPU else 'tnp_fitter_lowPU.py'
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--createBins'                   ])
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--createHists'                  ])
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--doFit'                        ])
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--doFit', '--mcSig',  '--altSig'])
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--doFit',             '--altSig'])
    cmds.append(['python3', ex, opt_e, opt_f, opt_m, '--sumUp'                        ])

    pretend = 0
    for cmd in cmds:
        if pretend:
            print(' '.join(cmd))
        else:
            subprocess.run(cmd)

mc    = 1
lowPU = 0

#flavs = ['el', 'mu'] if lowPU else ['mu']
flavs = ['mu', 'el'] if lowPU else ['mu']

for flav in flavs:

    #types = ['tracking', 'reco', 'idip', 'trigger', 'iso', 'isonotrig'] if flav == 'mu' else ['idiso', 'trigger']
    types = ['tracking']#, 'altreco', 'idip', 'iso', 'trigger', 'isonotrig'] if flav == 'mu' else ['idiso', 'trigger']

    for t in types:
    
        eras = ['BtoF', 'GtoH'] if not lowPU else ['lowPU'] 
        for e in eras :
    
            charges = ['plus', 'minus', 'both'] if t == 'trigger' else ['both']
        
            for c in charges:
    
                runCommands(flav,t,e,c,mc,lowPU)

