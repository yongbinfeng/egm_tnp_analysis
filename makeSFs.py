import ROOT, copy, datetime

today = datetime.date.today()

hists = []

prepost = {}

resdir = 'results_nodz_dxybs_mcTruth_2021-11-08_binnedInZ/'#nodz_dxybs_2021-09-14/' #'results_mcTruth' if mcTruth else 'results'
mcTruth = 'mcTruth' in resdir

for ch in ['plus', 'minus', 'both']:
    for era in ['lowPU', 'BtoF', 'GtoH', 'BtoH', 'B', 'C', 'D', 'E', 'F_preVFP', 'G', 'H']:
        for tnp in ['altreco', 'reco', 'tracking', 'idip', 'trigger', 'iso', 'isonotrig', 'alttrack', 'idiso']:
            for fl in ['mu', 'el']:
                tmp_infile = ROOT.TFile(resdir+'/efficiencies_{e}/{f}_{t}_{ch}/allEfficiencies_2D.root'.format(f=fl,e=era,ch=ch,t=tnp), 'read')
                if not tmp_infile.IsOpen():
                    continue
                print(tmp_infile)
                flstr = fl+'_' if era =='lowPU' else ''
                tmp_hist = tmp_infile.Get('SF2D_nominal')
                tmp_hist.SetName ('{f}SF2D_nominal_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_hist.SetTitle('{f}SF2D_nominal_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_da = tmp_infile.Get('EffData2D')
                tmp_da.SetName ('{f}effData_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_da.SetTitle('{f}effData_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_mc = tmp_infile.Get('EffMC2D')
                tmp_mc.SetName ('{f}effMC_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_mc.SetTitle('{f}effMC_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_daalteff = tmp_infile.Get('EffDataAltSig2D')
                tmp_daalteff.SetName ('{f}effData_altSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_daalteff.SetTitle('{f}effData_altSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_mcalteff = tmp_infile.Get('EffMCAltSig2D')
                tmp_mcalteff.SetName ('{f}effMC_altSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_mcalteff.SetTitle('{f}effMC_altSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_daalt = tmp_infile.Get('SF2D_dataAltSig')
                tmp_daalt.SetName ('{f}SF2D_dataAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_daalt.SetTitle('{f}SF2D_dataAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_mcalt = tmp_infile.Get('SF2D_MCAltSig')
                tmp_mcalt.SetName ('{f}SF2D_MCAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_mcalt.SetTitle('{f}SF2D_MCAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                tmp_damcalt = tmp_infile.Get('SF2D_dataMCAltSig')
                tmp_damcalt.SetName ('{f}SF2D_dataMCAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                tmp_damcalt.SetTitle('{f}SF2D_dataMCAltSig_{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                if era == 'BtoF':
                    prepost[f'da_{era}_{ch}_{tnp}'] = copy.deepcopy(tmp_da)
                    prepost[f'mc_{era}_{ch}_{tnp}'] = copy.deepcopy(tmp_mc)


                hists.extend([copy.deepcopy(tmp_hist), copy.deepcopy(tmp_da), copy.deepcopy(tmp_mc), copy.deepcopy(tmp_daalt), copy.deepcopy(tmp_mcalt), copy.deepcopy(tmp_damcalt), copy.deepcopy(tmp_daalteff), copy.deepcopy(tmp_mcalteff)])

                if tnp in ['iso', 'isonotrig']:

                    tmp_hist_data    = tmp_infile.Get('EffData2D')
                    tmp_hist_mc      = tmp_infile.Get('EffMC2D'  )
                    tmp_hist_altdata = tmp_infile.Get('EffDataAltSig2D')
                    tmp_hist_altmc   = tmp_infile.Get('EffMCAltSig2D'  )

                    tmp_hist_data    .SetName('{f}effData_anti{t}_{e}_{ch}'       .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_mc      .SetName('{f}effMC_anti{t}_{e}_{ch}'         .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_altdata .SetName('{f}effData_altSig_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_altmc   .SetName('{f}effMC_altSig_anti{t}_{e}_{ch}'  .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_data    .SetTitle('{f}effData_anti{t}_{e}_{ch}'       .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_mc      .SetTitle('{f}effMC_anti{t}_{e}_{ch}'         .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_altdata .SetTitle('{f}effData_altSig_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_hist_altmc   .SetTitle('{f}effMC_altSig_anti{t}_{e}_{ch}'  .format(e=era,ch=ch,t=tnp,f=flstr))

                    for ix in range(1,tmp_hist_data.GetNbinsX()+1):
                        for iy in range(1,tmp_hist_data.GetNbinsY()+1):
                            tmp_hist_data   .SetBinContent(ix, iy, 1.-tmp_hist_data   .GetBinContent(ix, iy))
                            tmp_hist_mc     .SetBinContent(ix, iy, 1.-tmp_hist_mc     .GetBinContent(ix, iy))
                            tmp_hist_altdata.SetBinContent(ix, iy, 1.-tmp_hist_altdata.GetBinContent(ix, iy))
                            tmp_hist_altmc  .SetBinContent(ix, iy, 1.-tmp_hist_altmc  .GetBinContent(ix, iy))

                    hists.append(copy.deepcopy(tmp_hist_data   ))
                    hists.append(copy.deepcopy(tmp_hist_mc     ))
                    hists.append(copy.deepcopy(tmp_hist_altdata))
                    hists.append(copy.deepcopy(tmp_hist_altmc  ))
                
                    tmp_sf_nominal  = tmp_hist_data   .Clone('{f}SF2D_nominal_anti{t}_{e}_{ch}'     .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_dataalt  = tmp_hist_altdata.Clone('{f}SF2D_dataAltSig_anti{t}_{e}_{ch}'  .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_mcalt    = tmp_hist_data   .Clone('{f}SF2D_MCAltSig_anti{t}_{e}_{ch}'    .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_datamcalt= tmp_hist_altdata.Clone('{f}SF2D_dataMCAltSig_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_nominal                 .SetTitle('{f}SF2D_nominal_anti{t}_{e}_{ch}'     .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_dataalt                 .SetTitle('{f}SF2D_dataAltSig_anti{t}_{e}_{ch}'  .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_mcalt                   .SetTitle('{f}SF2D_MCAltSig_anti{t}_{e}_{ch}'    .format(e=era,ch=ch,t=tnp,f=flstr))
                    tmp_sf_datamcalt               .SetTitle('{f}SF2D_dataMCAltSig_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp,f=flstr))

                    tmp_sf_nominal  .Divide(tmp_hist_mc)
                    tmp_sf_dataalt  .Divide(tmp_hist_mc)
                    tmp_sf_mcalt    .Divide(tmp_hist_altmc)
                    tmp_sf_datamcalt.Divide(tmp_hist_altmc)
                    #tmp_hist_data.SetName ('SF2D_nominal_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp))
                    #tmp_hist_data.SetTitle('SF2D_nominal_anti{t}_{e}_{ch}'.format(e=era,ch=ch,t=tnp))

                    hists.append(copy.deepcopy(tmp_sf_nominal  ))
                    hists.append(copy.deepcopy(tmp_sf_dataalt  ))
                    hists.append(copy.deepcopy(tmp_sf_mcalt    ))
                    hists.append(copy.deepcopy(tmp_sf_datamcalt))


                
                #outfile.Write(tmp_hist)

                #tmp_infile.Close()

of = 'allSFs_nodz_dxybs_mcTruth' if mcTruth else 'allSFs_nodz_dxybs'
if 'owPU' in resdir:
    of +='_lowPU'
outfile = ROOT.TFile(str(today)+'_'+of+'.root', 'recreate')
outfile.cd()

for h in hists:
    h.Write()

outfile.Close()
