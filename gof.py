import ROOT as rt
import os, optparse
import numpy as np
import uproot
import tdrstyle
tdrstyle.setTDRStyle()
rt.gROOT.SetBatch(True)
rt.gStyle.SetOptStat(0)
rt.gStyle.SetOptTitle(0)

lumi = 64
bin_edges = np.array([1200,1246,1313,1383,1455,1530,1607,1687,
                      1770,1856,1945,2037,2132,2231,2332,2438,
                      2546,2659,2775,2895,3019,3147,3279,3416,
                      3558,3704,3854,4010,4171,4337,4509,4686,
                      4869,5058,5253,5500,5663,5877,6099,6328,6564,6808]).astype('float') 

bin_edges = np.array([1313,1383,1455,1530,1607,1687,
                      1770,1856,1945,2037,2132,2231,2332,2438,
                      2546,2659,2775,2895,3019,3147,3279,3416,
                      3558,3704,3854,4010,4171,4337,4509,4686,
                      4869,5058,5253]).astype('float') 


n_bins = len(bin_edges)-1
max_bin = bin_edges[-1]
min_bin = bin_edges[0]
bws = (max_bin-min_bin)/n_bins

col = ['#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854','#ffd92f'] *3
col.reverse()

def load_data(indir,quantiles):

    histos_sig = {}
    histos_qcd = {}

    for q in quantiles:
        fname = os.path.join(indir, "data_mjj_sig_%s.root"%q)
        q_datafile = rt.TFile.Open(fname,'READ')
        tmp = q_datafile.Get("mjj_sig_%s"%q)
        histos_sig[q] = tmp.Rebin(n_bins,tmp.GetName()+"_dijetBins",bin_edges)
        histos_sig[q].SetDirectory(rt.gROOT)
        q_datafile.Close()

        fname = os.path.join(indir, "data_mjj_qcd_%s.root"%q)
        q_datafile = rt.TFile.Open(fname,'READ')
        tmp = q_datafile.Get("mjj_qcd_%s"%q)
        histos_qcd[q] = tmp.Rebin(n_bins,tmp.GetName()+"_dijetBins",bin_edges)
        histos_qcd[q].SetDirectory(rt.gROOT)
        q_datafile.Close()

    return histos_sig, histos_qcd

def make_ratio_plots(histos_sig, histos_qcd, quantiles, indir, mx, xs, width, legends): # Draw "nice" ratio plot for paper, check efficiency quantile/template and smoothness
    
    c_out = rt.TCanvas("ratio", "", 1000, 1000)
    pad1 = rt.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0.001)
    pad1.SetLeftMargin(0.13)
    pad1.Draw()
    pad1.cd()
    leg = rt.TLegend(0.52, 0.5, 0.89, 0.80)
    leg.SetNColumns(2)
    leg.SetTextSize(0.045)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    histos_qcd['q0'].SetLineColor(rt.kBlack)
    histos_qcd['q0'].SetMarkerColor(rt.kBlack)
    histos_qcd['q0'].Draw('HIST')
    histos_qcd['q0'].SetTitle("")
    histos_qcd['q0'].SetMinimum(0.2)
    histos_qcd['q0'].GetYaxis().SetTitleSize(0.05)
    histos_qcd['q0'].GetYaxis().SetLabelSize(0.05)
    histos_qcd['q0'].GetYaxis().SetTitleOffset(0.95)
    histos_qcd['q0'].SetMaximum(histos_qcd['q0'].Integral()*10.0)
    histos_qcd['q0'].GetYaxis().SetNdivisions(304)
    histos_sig['q90'].SetLineColor(rt.TColor.GetColor(col[-2]))
    histos_sig['q90'].SetFillColorAlpha((rt.TColor.GetColor(col[-2])), 0.30)
    histos_sig['q90'].Scale(10/histos_sig['q90'].Integral())
    histos_sig['q90'].Draw('HIST same')
    leg.SetHeader('Data:              G_{%s}(%.1f TeV, A.U.):'%(width,(mx/1000)))
    leg.AddEntry(histos_qcd['q0']  ,legends[0], 'lep')
    leg.AddEntry(histos_sig['q90']  ,legends[-2], 'l')

    for i,q in enumerate(quantiles):
        if q == 'q0':
            continue
        histos_qcd[q].SetLineColor(rt.TColor.GetColor(col[i]))
        histos_qcd[q].SetMarkerColor(rt.TColor.GetColor(col[i]))
        histos_qcd[q].Draw('pez same')
        leg.AddEntry(histos_qcd[q]  , legends[i], 'lep')
        leg.AddEntry(None, '', '')
    histos_qcd['q0'].Draw('HISTsame')
    leg.Draw('same')
    pad1.SetLogy()
    pad1.SetTitle('')
    c_out.cd()
    pad2 = rt.TPad("pad2", "pad2", 0, 0, 1, 0.3)
    pad2.SetTopMargin(0.00)
    pad2.SetBottomMargin(0.27)
    pad2.SetLeftMargin(0.13)
    # pad2.SetGrid()
    pad2.Draw()
    pad2.cd()
    data_hist_ratio = {}
    tline = {}
    for i,q in enumerate(quantiles):
        histos_qcd[q].Sumw2()
        data_hist_ratio[q] = histos_qcd[q].Clone('data_hist_ratio_{}'.format(i))
        data_hist_ratio[q].Divide(histos_qcd['q0'])
        data_hist_ratio[q].Scale(histos_qcd['q0'].Integral()/histos_qcd[q].Integral())
        print("Ratio for quantile {Q} is {R}".format(Q=q,R=histos_qcd['q0'].Integral()/histos_qcd[q].Integral()))
        data_hist_ratio[q].binning = bin_edges
        if i == 1:
            data_hist_ratio[q].SetTitle('')
            data_hist_ratio[q].Draw('pez')
            data_hist_ratio[q].SetMaximum(1.6)
            data_hist_ratio[q].SetMinimum(0.4)
            data_hist_ratio[q].SetYTitle('#frac{#epsilon#timesq}{q_{0.0-0.3}}')
            data_hist_ratio[q].GetYaxis().CenterTitle()
            data_hist_ratio[q].SetXTitle('M_{jj} (GeV)')
            data_hist_ratio[q].GetYaxis().SetTitleOffset(0.4)
            data_hist_ratio[q].GetYaxis().SetTitleSize(0.12)
            data_hist_ratio[q].GetYaxis().SetLabelSize(0.12)
            data_hist_ratio[q].GetYaxis().SetNdivisions(304)
            data_hist_ratio[q].GetXaxis().SetNdivisions(909)
            data_hist_ratio[q].GetXaxis().SetTitleOffset(0.95)
            data_hist_ratio[q].GetXaxis().SetTitleSize(0.12)
            data_hist_ratio[q].GetXaxis().SetLabelSize(0.12)
            data_hist_ratio[q].GetXaxis().SetTickSize(0.07)
        else:
            data_hist_ratio[q].Draw('pez same')
        tline[q] = rt.TLine(min_bin, 1.0, max_bin, 1.0)
        tline[q].SetLineColor(rt.kBlack)
        tline[q].SetLineStyle(rt.kDashed)
        tline[q].Draw('same')

    c_out.cd()
    latex = rt.TLatex()
    latex.SetNDC ()
    latex.SetTextSize (0.033)
    latex.SetTextFont( 42 )
    latex.DrawLatex (0.16 ,0.88 , "%.0f fb^{-1} dijet events"%lumi)
    c_out.Update()
    c_out.Draw()
    c_out.SaveAs(os.path.join(indir, 'mjj_gof_ratio.pdf'))
    
def makeWS(data_rej, data_acc, signal_rej, signal_acc, outdir, quantile):

    basedir = os.getcwd()
    workdir = '{}/{}'.format(basedir,outdir)
    os.chdir(workdir)

    outname  = 'datacard_gof_%s.root'%quantile
    efficiency = data_acc.Integral()/data_rej.Integral() #How much rej histogram needs to be scaled to match acc
    datacard_ws = rt.TFile.Open(outname,'recreate')
    w = rt.RooWorkspace('w','w')
    x = rt.RooRealVar('x','x',min_bin,max_bin)
    w.factory('x[%.1f,%.1f]'%(min_bin, max_bin))
    acc_bin_functions = rt.RooArgList()
    rej_bin_functions = rt.RooArgList()
    
    w.factory('eff_%s[%f,%f,%f]'%(quantile,efficiency,efficiency*0.95,efficiency*1.05)) # Name efficiency per quantile to allow for combination
    w.var('eff_%s'%quantile).setConstant(False)
    
    empty_hist = rt.TH1D('empty_hist','empty_hist', n_bins, bin_edges.astype('float'))

    for iBinX in range(1,n_bins+1):
        empty_hist.SetBinContent(iBinX,1)
        rej_bin = data_rej.GetBinContent(iBinX)
        w.factory('crBin%i_In[%.1f]'%(iBinX,rej_bin))
        w.var('crBin%i_In'%iBinX).setConstant(True) # Fix control region bin count
        w.factory('crBin_%i[0,-10,10]'%(iBinX)) # How much data is pulled away from the background model
        w.var('crBin_%i'%(iBinX)).setConstant(False)
        if rej_bin !=  0.:
            power = 1/rt.TMath.Sqrt(rej_bin)
        else: # Skip zero bins
            power = -1.0
            w.var('crBin_%i'%(iBinX)).setConstant(True)
            
        w.factory("expr::crBin%iFunc('max(0,@0*pow(1.0+%f,@1))',crBin%i_In,crBin_%i)"%(iBinX,power,iBinX,iBinX)) 
        w.factory("expr::bin%iFunc_q%s('max(0,@0*@1)',eff_%s,crBin%iFunc)"%(iBinX,quantile,quantile,iBinX))
        rej_bin_functions.add(w.function('crBin%iFunc'%(iBinX)))
        acc_bin_functions.add(w.function('bin%iFunc_q%s'%(iBinX,quantile)))
        
    qcd_rph_rej = rt.RooParametricHist('background_rej','background_rej',w.var('x'),rej_bin_functions,empty_hist)
    qcd_rph_rej_norm = rt.RooAddition('background_rej_norm','background_rej_norm',rej_bin_functions)
    qcd_rph_acc = rt.RooParametricHist('background_acc','background_acc',w.var('x'),acc_bin_functions,empty_hist)
    qcd_rph_acc_norm = rt.RooAddition('background_acc_norm','background_acc_norm',acc_bin_functions)
    getattr(w,'import')(qcd_rph_rej, rt.RooCmdArg())
    getattr(w,'import')(qcd_rph_rej_norm, rt.RooFit.RecycleConflictNodes())
    getattr(w,'import')(qcd_rph_acc, rt.RooCmdArg())
    getattr(w,'import')(qcd_rph_acc_norm, rt.RooFit.RecycleConflictNodes())

    ds_signal_acc = rt.RooDataHist('signal_acc','signal_acc',rt.RooArgList(w.var('x')),signal_acc)
    ds_signal_rej = rt.RooDataHist('signal_rej','signal_rej',rt.RooArgList(w.var('x')),signal_rej)
    getattr(w,'import')(ds_signal_acc, rt.RooCmdArg())
    getattr(w,'import')(ds_signal_rej, rt.RooCmdArg())

    ds_data_acc = rt.RooDataHist('data_obs_acc','data_obs_acc',rt.RooArgList(w.var('x')),data_acc) #
    ds_data_rej = rt.RooDataHist('data_obs_rej','data_obs_rej',rt.RooArgList(w.var('x')),data_rej)
    getattr(w,'import')(ds_data_acc, rt.RooCmdArg())
    getattr(w,'import')(ds_data_rej, rt.RooCmdArg())

    datacard_ws.cd()
    w.Write()
    datacard_ws.Close()

    datacard_ratio = \
    '''
    imax 1
    jmax 1
    kmax *
    ---------------
    shapes * * {WS} w:$PROCESS_$CHANNEL w:$PROCESS_$CHANNEL_$SYSTEMATIC
    ---------------
    bin {BIN}
    observation {OBS}
    ------------------------------
    bin             {BIN}      {BIN}
    process         signal     background
    process         0          1
    rate            {SIGRATE}    {BKGRATE}
    --------------------------------
    lumi lnN 1.01 -
    '''
    datacard_ratio += 'eff_%s   flatParam\n'%quantile

    for i in range(1,n_bins+1):
        datacard_ratio += 'crBin_%i   flatParam\n'%(i)

    datacard_ratio_acc = datacard_ratio.format(BIN='acc',
                                OBS=-1,
                                BKGRATE=1,
                                SIGRATE=signal_acc.Integral(),
                                WS=outname)
    with open(outname.replace('.root','_acc.txt'),'w') as f:
        f.write(datacard_ratio_acc)
        
        
    datacard_ratio_rej = datacard_ratio.format(BIN='rej',
                                OBS=-1,
                                BKGRATE=1,
                                SIGRATE=signal_rej.Integral(),
                                WS=outname)
    with open(outname.replace('.root','_rej.txt'),'w') as f:
        f.write(datacard_ratio_rej)
    os.system('combineCards.py rej={REJ} acc={ACC} > {RATIO}'.format(REJ=outname.replace('.root','_rej.txt'),ACC=outname.replace('.root','_acc.txt'),RATIO=outname.replace('.root','_ratio.txt')))
    os.chdir(basedir)

def plotGOF(obs_gof, exp_gof, quantile, n_dof=n_bins):
    
    n_extreme = len(exp_gof[exp_gof > obs_gof])
    n_total = len(exp_gof)
    pval_toys = 1.*n_extreme/n_total
    pval = rt.TMath.Prob(obs_gof,n_dof)  # get p-value assuming chi2 dist (may not be valid)
    
    bin_width = (max(exp_gof+[obs_gof])+np.std(exp_gof)-(min(exp_gof)-np.std(exp_gof)))/30.
    exp_gof_hist = rt.TH1D('gof','gof',30,min(exp_gof)-np.std(exp_gof), max(exp_gof+[obs_gof])+np.std(exp_gof))
    exp_gof_hist_gt = rt.TH1D('gof_gt','gof_gt',30,min(exp_gof)-np.std(exp_gof), max(exp_gof+[obs_gof])+np.std(exp_gof))
    for g in exp_gof:
        exp_gof_hist.Fill(g)
        if g > obs_gof:
            exp_gof_hist_gt.Fill(g)
    
    d = rt.TCanvas("ratio", "", 1000, 800)
    d.SetLeftMargin(0.13)
    f = rt.TF1("chi2","%f*ROOT::Math::chisquared_pdf(x,%i,0)"%(exp_gof_hist.Integral()*bin_width,n_dof),min(exp_gof)-np.std(exp_gof),max(exp_gof+[obs_gof])+np.std(exp_gof))
    tleg = rt.TLegend(0.48, 0.6, 0.89, 0.85)
    exp_gof_hist.Draw('hist')
    exp_gof_hist.SetXTitle('Test statistic -2ln#lambda')
    exp_gof_hist.SetYTitle('N toys')
    exp_gof_hist.SetTitle("")
    exp_gof_hist.GetYaxis().SetLabelSize(0.05)
    exp_gof_hist.GetYaxis().SetTitleSize(0.05)
    f.SetLineColor((rt.TColor.GetColor(col[1])))
    exp_gof_hist.SetLineWidth(2)
    exp_gof_hist.SetLineColor((rt.TColor.GetColor(col[0])))
    exp_gof_hist_gt.SetLineColor((rt.TColor.GetColor(col[0])))
    exp_gof_hist_gt.SetFillColorAlpha((rt.TColor.GetColor(col[0])), 0.30)
    exp_gof_hist_gt.Draw('fhistsame')
    f.Draw("same")
    line = rt.TLine(obs_gof,0,obs_gof,exp_gof_hist.GetMaximum())
    line.SetLineWidth(2)
    line.Draw()
    tleg.AddEntry(exp_gof_hist_gt,'p-value (from toys) = %.4f'%pval_toys)
    tleg.AddEntry(f,'p-value (from #chi^{2}) = %.4f'%pval,'l')
    tleg.AddEntry(line,'Observed test stat. = {:.1f}'.format(obs_gof),'l')
    tleg.Draw()
    d.Draw()
    d.Update()
    d.SaveAs("GOF_{q}.pdf".format(q=quantile))

def runFitDiagnosis(datacard, quantile, cats=['rej', 'acc']):
    
    os.system('combine -M FitDiagnostics -d {DATACARD} -n _{Q} --saveShapes --saveWithUncertainties --dataset data_obs --verbose 0'.format(DATACARD=datacard, Q=quantile))
    fitDiag = rt.TFile.Open('fitDiagnostics_{Q}.root'.format(Q=quantile),'r')

    f = rt.TCanvas('f','f',1000,800)
    f.cd()

    tleg = rt.TLegend(0.48, 0.6, 0.89, 0.85)
    tleg.SetTextSize(0.05)
    tleg.SetBorderSize(0)
    tleg.SetFillStyle(0)
    tleg.SetTextSize (0.03)
    tleg.SetTextFont( 62 )
    tleg.SetTextSize (0.03)
    tleg.SetTextFont( 42 )

    byhand_gof = 0.0

    for cat in cats:
        print cat
        bkgd = fitDiag.Get('shapes_fit_b/{cat}/background'.format(cat=cat))
        bkgd.Scale(bws) # need to multiply by bin width for some reason?
        data = fitDiag.Get('shapes_fit_b/{cat}/data'.format(cat=cat))
        if cat.find('rej')!=-1:
            bkgd.Draw('hist')
            bkgd.SetLineColor(rt.kRed)
            tleg.AddEntry(bkgd,'rej','l')
        else:
            bkgd.Draw("histsame")
            tleg.AddEntry(bkgd,'acc','l')
        data.SetMarkerStyle(20)
        data.SetMarkerColor(rt.kBlack)
        for i in range(0,bkgd.GetNbinsX()):
            data.SetPointEXlow(i,0)
            data.SetPointEXhigh(i,0)
            data.SetPoint(i,data.GetX()[i], bws*data.GetY()[i]) # need to multiply by bin width for some reason?
            data.SetPointEYlow(i,bws*data.GetErrorYlow(i)) # need to multiply by bin width for some reason?
            data.SetPointEYhigh(i,bws*data.GetErrorYhigh(i)) # need to multiply by bin width for some reason?
        data.Draw('samepez')

        for i in range(0,bkgd.GetNbinsX()):
            mjjvalue = bkgd.GetBinCenter(i+1)
            fi = bkgd.GetBinContent(i+1)
            di = data.GetY()[i]
            if di == 0:
                print ("No data at mjj", mjjvalue)
                print ("fi", fi)
                print ("di", di)
                gofi = 0.
            else:    
                gofi = 2*(fi - di + di*rt.TMath.Log(di/fi)) # see eq. 14 of http://cousins.web.cern.ch/cousins/ongoodness6march2016.pdf # expect each bin to give GOF contribution ~ O(1)
            if gofi>5:
                print(" Category {}: mjj bin = {} gof = {} (data = {} fit = {})".format(cat,mjjvalue,gofi, di, fi))
            byhand_gof += gofi

    tleg.Draw('same')
    f.SetLogy()
    f.Draw()
    f.SaveAs("fitResults_{}.pdf".format(quantile))
    print("By hand obs GOF = {}".format(byhand_gof))

def runCombine(quantile,datacarddir):

    basedir = os.getcwd()
    workdir = '{}/{}'.format(basedir,datacarddir)
    cardname  = os.path.join(workdir, 'datacard_gof_%s_ratio.txt'%quantile)

    os.chdir(workdir)
    os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d {DATACARD} -n gof_{Q} --dataset data_obs -v 0'.format(DATACARD=cardname, Q=quantile))
    
    if options.runToys:

        toys_per_job = int(options.ntoys)/2
        print("Running %s toys with 5 different seeds!"%toys_per_job)
        for i in range(2):
            os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d {DATACARD}  -t {NTOYS} --toysFreq -n gof_toys_{Q}  --dataset data_obs -s {S} -v 0'.format(DATACARD=cardname, Q=quantile, NTOYS=toys_per_job,S=40+i))
        os.system('hadd -f higgsCombinegof_toys_{Q}.GoodnessOfFit.mH120.ALLTOYS.root higgsCombinegof_toys_{Q}.GoodnessOfFit.mH120.4*.root'.format(Q=quantile))
    
        obs_gof_file = uproot.open('higgsCombinegof_{Q}.GoodnessOfFit.mH120.root'.format(Q=quantile))
        obs_gof = obs_gof_file['limit'].arrays('limit')['limit'][0]
        exp_gof_file = uproot.open('higgsCombinegof_toys_{Q}.GoodnessOfFit.mH120.ALLTOYS.root'.format(Q=quantile))
        exp_gof = exp_gof_file['limit'].arrays('limit')['limit']
        print("Obs.   {:.1f}".format(obs_gof))
        print("Exp.   {:.1f}\n".format(np.mean(exp_gof)))   
        runFitDiagnosis(cardname, quantile)
        plotGOF(obs_gof,exp_gof,q)

    os.chdir(basedir)

def runCombination(datacarddir, qacc=['q30', 'q50', 'q70', 'q90']):
  
  print "\n Doing full combination of all quantiles"

  basedir = os.getcwd()
  workdir = '{}/{}'.format(basedir,datacarddir)
  card_prefix  = os.path.join(workdir, 'datacard_gof_')
  
  combined_card_name = 'combined_gof.txt'
  os.chdir(workdir)
  os.system('rm {CCARD}'.format(CCARD=combined_card_name))

  command = 'combineCards.py rej={CARD}{Q}_rej.txt'.format(CARD=card_prefix,Q=qacc[0])
 
  for q in qacc:
    command += ' {Q}={CARD}{Q}_acc.txt'.format(CARD=card_prefix, Q=q)
  command += ' &> {CCARD}'.format(CCARD=combined_card_name)
  print("Combining:", command)

  os.system(command)
  os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d {CCARD} -n gof_combined --dataset data_obs -v 3'.format(CCARD=combined_card_name))
  
  if options.runToys:
    toys_per_job = int(options.ntoys)/5
    print("Running %s toys with 5 different seeds!"%toys_per_job)
    for i in range(5):
        os.system('combine -M GoodnessOfFit --algo saturated --fixedSignalStrength 0 -d {CCARD}  -t {NTOYS} --toysFreq -n gof_toys_combined --dataset data_obs -v 0 -s {S}'.format(CCARD=combined_card_name,NTOYS=toys_per_job,S=22+i))
    os.system('hadd -f higgsCombinegof_toys_combined.GoodnessOfFit.mH120.ALLTOYS.root higgsCombinegof_toys_combined.GoodnessOfFit.mH120.4*.root')
    obs_gof_file = uproot.open('higgsCombinegof_combined.GoodnessOfFit.mH120.root')
    obs_gof = obs_gof_file['limit'].arrays('limit')['limit'][0]
    exp_gof_file = uproot.open('higgsCombinegof_toys_combined.GoodnessOfFit.mH120.ALLTOYS.root')
    exp_gof = exp_gof_file['limit'].arrays('limit')['limit']
    print("Obs.   {:.1f}".format(obs_gof))
    print("Exp.   {:.1f}\n".format(np.mean(exp_gof)))   
    runFitDiagnosis(combined_card_name, "COMBINED", cats=['rej','q30', 'q50', 'q70', 'q90'])
    plotGOF(obs_gof,exp_gof,"COMBINED", n_dof=n_bins*len(qacc))
  os.chdir(basedir)
 


if __name__ == "__main__":

   #python gof.py -i inputdir

   #some configuration
   parser = optparse.OptionParser()
   parser.add_option("--xsec","--xsec",dest="xsec",type=float,default=0.0006,help="Injected signal cross section in fb")
   parser.add_option("-M","-M",dest="mass",type=float,default=3500.,help="Injected signal mass")
   parser.add_option("-i","--indir",dest="indir",default='./',help="directory with histogram outputs from dijet fitting code")
   parser.add_option("--res", "--res", dest="sig_res", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
   parser.add_option('-N','--ntoys', dest='ntoys',default=100)
   parser.add_option('-T','--runToys',action='store_true', dest='runToys'  , default=False, help='runToys')
   parser.add_option('-C','--doCombination',action='store_true', dest='doCombination'  , default=False, help='do full combination')
   parser.add_option('-Q','--perQuantile',action='store_true', dest='doPerQuantile'  , default=False, help='do per quantile fit')

   (options,args) = parser.parse_args()
    
   mass = options.mass
   sig_res = options.sig_res
   indir = options.indir
   outdir = indir

   quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'total'] #Delphes inverted wrt CASE -- most anomalous is q90
   legends  = ['q = 0-30', 'q = 30-50', 'q = 50-70', 'q = 70-90', 'q = 90-100', 'Inclusive'] #Delphes inverted wrt CASE -- most anomalous is q90

   histos_sig, histos_qcd = load_data(indir,quantiles)
   make_ratio_plots(histos_sig, histos_qcd, quantiles, indir, mass, options.xsec, options.sig_res, legends)
   
   for q in quantiles:
    if q not in ['q0','total']:
        print("Making workspaces for {}".format(q))
        # makeWS(histos_qcd['q0'], histos_qcd[q], histos_sig['q0'], histos_sig[q], outdir, q)
        # runCombine(q,outdir)

   if options.doCombination:
      runCombination(datacarddir=outdir, qacc=['q30', 'q50', 'q70', 'q90'])