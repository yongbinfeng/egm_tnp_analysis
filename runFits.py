import subprocess
#from multiprocessing import Process
import time
import argparse

def runCommands(wp, era, inputMC, inputData, options):
    outdir = options.outdir
    pretend = options.dryRun
    print()
    print("-"*30)
    print("Working point = ",wp)
    print("-"*30)
    opt_e = '--era='+era
    opt_f = '--flag='+wp
    opt_iMC = '--inputMC='+inputMC
    opt_iData = '--inputData='+inputData
    cmds = []
    ex = 'tnpEGM_fitter.py'
    cmds.append(['python', ex, opt_e, opt_f, '--createBins'                   ])
    cmds.append(['python', ex, opt_e, opt_f, opt_iMC , opt_iData , '--createHists'])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit',                        ])
    #cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig'                        ])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig',  '--altSig'])
    cmds.append(['python', ex, opt_e, opt_f, '--doFit',             '--altSig'])
    #cmds.append(['python', ex, opt_e, opt_f, '--doFit', '--mcSig',  '--altBkg'])
    #cmds.append(['python', ex, opt_e, opt_f, '--doFit',             '--altBkg'])
    cmds.append(['python', ex, opt_e, opt_f, '--sumUp'                        ])

    for cmd in cmds:
        if outdir:
            cmd.append(f"--outdir={outdir}")
        if options.useTrackerMuons:
            cmd.append("--useTrackerMuons")
        if pretend:
            print(' '.join(cmd))
        else:
            subprocess.run(cmd, check=True)

working_points_global = {
    ## for global muons
    ##
    'mu_reco_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_reco_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_reco_data_vertexWeights1_oscharge1.root'],
    # 'mu_reco_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_recoplus_mc_vertexWeights1_oscharge1.root',
    #                  '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_recoplus_data_vertexWeights1_oscharge1.root'],
    # 'mu_reco_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_recominus_mc_vertexWeights1_oscharge1.root',
    #                  '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_recominus_data_vertexWeights1_oscharge1.root'],
    # OK
    'mu_tracking_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_tracking_pt25to65_mass50to130_4ptbins/tnp_tracking_mc_vertexWeights1_oscharge0.root',
                         '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_tracking_pt25to65_mass50to130_4ptbins/tnp_tracking_data_vertexWeights1_oscharge0.root'],
    'mu_idip_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_idip_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_idip_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_triggerplus_mc_vertexWeights1_oscharge1.root',
                        '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_triggerplus_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_triggerminus_mc_vertexWeights1_oscharge1.root',
                        '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_triggerminus_data_vertexWeights1_oscharge1.root'],
    #'mu_iso_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_iso_mc_vertexWeights1_oscharge1.root',
    #                '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_iso_data_vertexWeights1_oscharge1.root'],
    'mu_iso_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_byCharge/tnp_isoplus_mc_vertexWeights1_oscharge1.root',
                    '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_byCharge/tnp_isoplus_data_vertexWeights1_oscharge1.root'],
    'mu_iso_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_byCharge/tnp_isominus_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/test_globalMuons_byCharge/tnp_isominus_data_vertexWeights1_oscharge1.root'],
    'mu_isonotrig_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_isonotrig_mc_vertexWeights1_oscharge1.root',
                          '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_isonotrig_data_vertexWeights1_oscharge1.root'],
    'mu_veto_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_veto_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/newHisto_allWP_mass60to120_noTracking/tnp_veto_data_vertexWeights1_oscharge1.root'],
}

working_points_tracker = {
    'mu_reco_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_reco_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_reco_data_vertexWeights1_oscharge1.root'],
    'mu_tracking_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_tracking_mc_vertexWeights1_oscharge0.root',
                         '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_tracking_data_vertexWeights1_oscharge0.root'],
    'mu_idip_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_idip_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_idip_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_plus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_triggerplus_mc_vertexWeights1_oscharge1.root',
                        '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_triggerplus_data_vertexWeights1_oscharge1.root'],
    'mu_trigger_minus': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_triggerminus_mc_vertexWeights1_oscharge1.root',
                        '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_triggerminus_data_vertexWeights1_oscharge1.root'],
    'mu_iso_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_iso_mc_vertexWeights1_oscharge1.root',
                    '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_iso_data_vertexWeights1_oscharge1.root'],
    'mu_isonotrig_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_isonotrig_mc_vertexWeights1_oscharge1.root',
                          '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_isonotrig_data_vertexWeights1_oscharge1.root'],
    'mu_veto_both': ['/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_veto_mc_vertexWeights1_oscharge1.root',
                     '/home/m/mciprian/tnp/Steve_Marc_Raj/outputs/trackerMuons_highPurity_allWP/tnp_veto_data_vertexWeights1_oscharge1.root'],
}


#working_points = {'mu_iso_both': ['/home/users/rajarshi/Steve_Erc/Isolation_MC_full_2_08_2022.root','/home/users/rajarshi/Steve_Erc/Isolation_Data_full_2_08_2022.root'],}
#for i in range(1,41):
    #working_points['mu_trigger_both_qtbin{0}'.format(i)] = ['input_MC_file_{0}.root'.format(i),'input_Data_file_{0}.root'.format(i)] 
#    working_points['mu_iso_both_qtbin{0}'.format(i)] = ['/home/users/rajarshi/Steve_Erc/Isolation_MC_2_08_2022_{0}.root'.format(i),'/home/users/rajarshi/Steve_Erc/Isolation_Data_2_08_2022_{0}.root'.format(i)] 


parser = argparse.ArgumentParser()
parser.add_argument('-o',  '--outdir', default=None, type=str,
                    help='name of the output folder (if not passed, a default one is used, which has the time stamp in it)')
parser.add_argument('-e',  '--era', default=['GtoH'], nargs='+', type=str, choices=['GtoH', 'BtoF'],
                    help='Choose the era')
parser.add_argument('-d',  '--dryRun', action='store_true',
                    help='Do not execute commands, just print them')
parser.add_argument('-s','--steps', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in working_points_global.keys()]),
                    help='Default runs all working points, but can choose only some if needed')
parser.add_argument('-x','--exclude', default=None, nargs='*', type=str, choices=list([x.split("_")[1] for x in working_points_global.keys()]),
                    help='Default runs all working points, but can exclude some if needed')
parser.add_argument('--useTrackerMuons', action='store_true'  , help = 'Measuring efficiencies specific for tracker muons (different tunings needed')
args = parser.parse_args()

tstart = time.time()
cpustrat = time.process_time()

if args.exclude and args.steps:
    print("Warning: --exclude and --steps are not supposed to be used together. Try again")
    quit()

eras = args.era
working_points = working_points_tracker if args.useTrackerMuons else working_points_global

stepsToRun = []
if args.steps:
    for x in working_points.keys():
        step = x.split("_")[1]
        if step in args.steps:
            stepsToRun.append(x)
else:
    if args.exclude:
        for x in working_points.keys():
            step = x.split("_")[1]
            if step not in args.exclude:
                stepsToRun.append(x)
    else:
        stepsToRun = working_points.keys()        
    
#procs = []
for e in eras:
    for wp in working_points:
        if wp not in stepsToRun:
            continue
        inputMC = working_points[wp][0]
        inputData = working_points[wp][1]
        runCommands( wp, e, inputMC, inputData, args)

elapsed = time.time() - tstart
elapsed_cpu = time.process_time() - cpustrat
print()
print()
print('Execution time:', elapsed, 'seconds')
print('CPU Execution time:', elapsed_cpu , 'seconds')
print()
