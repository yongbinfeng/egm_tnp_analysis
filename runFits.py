import subprocess
#from multiprocessing import Process

def runCommands(wp,era,inputMC,inputData,mc=False,lowPU=False):
    print("Working point = ",wp)
    opt_e = '--era='+era
    opt_f = '--flag='+wp
    #opt_m = '--mcTruth' if mc else ''
    opt_iMC = '--inputMC='+inputMC
    opt_iData = '--inputData='+inputData
    cmds = []
    ex = 'tnpEGM_fitter.py' if not lowPU else 'tnp_fitter_lowPU.py'
    cmds.append(['python', ex, opt_e, opt_f, '--createBins'                   ])
    cmds.append(['python', ex, opt_e, opt_f, opt_iMC , opt_iData , '--createHists'])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit'                        ])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig'                        ])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig',  '--altSig'])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit',             '--altSig'])
    cmds.append(['python', ex, opt_e, opt_f, '--sumUp'                        ])

    pretend = 0
    for cmd in cmds:
        if pretend:
            print(' '.join(cmd))
        else:
            subprocess.run(cmd)

#mc    = 1
#lowPU = 0

#flavs = ['el', 'mu'] if lowPU else ['mu']
#flavs = ['mu', 'el'] if lowPU else ['mu']

#for flav in flavs:

#    types = ['tracking', 'reco', 'idip', 'trigger', 'iso'] if flav == 'mu' else ['idiso', 'trigger']
#    types = ['tracking']#, 'altreco', 'idip', 'iso', 'trigger', 'isonotrig'] if flav == 'mu' else ['idiso', 'trigger']

#    for t in types:
    
#        eras = ['BtoF', 'GtoH'] if not lowPU else ['lowPU'] 
#        for e in eras :
    
#            charges = ['both_qtbin{0}'.format(i) for i in range(1,41)] if (t == 'trigger' or t == 'iso') else ['both']
        
#            for c in charges:

#working_points = {
#     'mu_reco_both': ['input_MC_file','input_Data_file'],
#     'mu_tracking_both': ['input_MC_file','input_Data_file'],
#     'mu_idip_both': ['input_MC_file','input_Data_file'],
#     'mu_trigger_both': ['input_MC_file','input_Data_file'],
#     'mu_iso_both': ['input_MC_file','input_Data_file']
#}

import time

tstart = time.time()
cpustrat = time.process_time()


working_points = {'mu_iso_both': ['/home/users/rajarshi/Steve_Erc/Isolation_MC_full_2_08_2022.root','/home/users/rajarshi/Steve_Erc/Isolation_Data_full_2_08_2022.root'],}
for i in range(1,41):
    #working_points['mu_trigger_both_qtbin{0}'.format(i)] = ['input_MC_file_{0}.root'.format(i),'input_Data_file_{0}.root'.format(i)] 
    working_points['mu_iso_both_qtbin{0}'.format(i)] = ['/home/users/rajarshi/Steve_Erc/Isolation_MC_2_08_2022_{0}.root'.format(i),'/home/users/rajarshi/Steve_Erc/Isolation_Data_2_08_2022_{0}.root'.format(i)] 
     
eras = ['GtoH']

#procs = []
for e in eras:
    for wp in working_points:
        inputMC = working_points[wp][0]
        inputData = working_points[wp][1]
        runCommands(wp,e,inputMC,inputData)
        #proc = Process(target=runCommands, args=(wp,e,inputMC,inputData,))
        #procs.append(proc)
        #proc.start()
        

#for proc in procs:
#    proc.join()


elapsed = time.time() - tstart
elapsed_cpu = time.process_time() - cpustrat
print('Execution time:', elapsed, 'seconds')
print('CPU Execution time:', elapsed_cpu , 'seconds')
