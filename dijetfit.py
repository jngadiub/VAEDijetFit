import h5py, math, commands, random
from array import array
import numpy as np
import time, sys, os, optparse, json

import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.RooRandom.randomGenerator().SetSeed(random.randint(0, 1e+6))

from Fitter import Fitter
from DataCardMaker import DataCardMaker
from Utils import *
 
def get_generated_events(filename):

 with open('files_count.json') as f:
  data = json.load(f)
 
 N = 0
 for k in data.keys():
  if (k in filename or k.replace('_EXT','') in filename) and not 'SIDEBAND' in k: N+=data[k][0]
  
 return N 
                                                               
def PlotFitResults(frame,fitErrs,nPars,pulls,data_name,pdf_name,chi2,ndof,canvname):

 c1 =ROOT.TCanvas("c1","",800,800)
 c1.SetLogy()
 c1.Divide(1,2,0,0,0)
 c1.SetLogy()
 c1.cd(1)
 p11_1 = c1.GetPad(1)
 p11_1.SetPad(0.01,0.26,0.99,0.98)
 p11_1.SetLogy()
 p11_1.SetRightMargin(0.05)

 p11_1.SetTopMargin(0.1)
 p11_1.SetBottomMargin(0.02)
 p11_1.SetFillColor(0)
 p11_1.SetBorderMode(0)
 p11_1.SetFrameFillStyle(0)
 p11_1.SetFrameBorderMode(0)
 frame.GetYaxis().SetTitleSize(0.06)
 frame.GetYaxis().SetTitleOffset(0.98)
 frame.SetMinimum(0.2)
 frame.SetMaximum(1E7)
 frame.SetName("mjjFit")
 frame.GetYaxis().SetTitle("Events / 100 GeV")
 frame.SetTitle("")
 frame.Draw()
    
 legend = ROOT.TLegend(0.45097293,0.64183362,0.6681766,0.879833)
 legend2 = ROOT.TLegend(0.45097293,0.64183362,0.6681766,0.879833)
 legend.SetTextSize(0.046)
 legend.SetLineColor(0)
 legend.SetShadowColor(0)
 legend.SetLineStyle(1)
 legend.SetLineWidth(1)
 legend.SetFillColor(0)
 legend.SetFillStyle(0)
 legend.SetMargin(0.35)
 legend2.SetTextSize(0.038)
 legend2.SetLineColor(0)
 legend2.SetShadowColor(0)
 legend2.SetLineStyle(1)
 legend2.SetLineWidth(1)
 legend2.SetFillColor(0)
 legend2.SetFillStyle(0)
 legend2.SetMargin(0.35)
 legend.AddEntry(frame.findObject(data_name),"Data","lpe")
 legend.AddEntry(frame.findObject(pdf_name),"%i par. background fit"%nPars,"l")
 legend2.AddEntry("","","")
 legend2.AddEntry("","","")
 legend2.AddEntry("","","")
 legend2.AddEntry("","","")
 legend2.AddEntry(frame.findObject(fitErrs),"","f")
 legend2.AddEntry("","","")

 legend2.Draw("same")
 legend.Draw("same")

 pt = ROOT.TPaveText(0.18,0.06,0.54,0.17,"NDC")
 pt.SetTextFont(42)
 pt.SetTextAlign(12)
 pt.SetFillColor(0)
 pt.SetBorderSize(0)
 pt.SetFillStyle(0)
 pt.AddText("Chi2/ndf = %.2f/%i = %.2f"%(chi2,ndof,chi2/ndof))
 pt.AddText("Prob = %.3f"%ROOT.TMath.Prob(chi2,ndof))
 pt.Draw()
 
 c1.Update()

 c1.cd(2)
 p11_2 = c1.GetPad(2)
 p11_2.SetPad(0.01,0.02,0.99,0.27)
 p11_2.SetBottomMargin(0.35)
 p11_2.SetRightMargin(0.05)
 p11_2.SetGridx(0)
 p11_2.SetGridy(0)
 pulls.SetMinimum(-10)
 pulls.SetMaximum(10)
 pulls.SetTitle("")
 pulls.SetXTitle("Dijet invariant mass (GeV)")
 pulls.GetXaxis().SetTitleSize(0.06)
 pulls.SetYTitle("#frac{Data-Fit}{#sigma_{data}}")
 pulls.GetYaxis().SetTitleSize(0.15)
 pulls.GetYaxis().CenterTitle()
 pulls.GetYaxis().SetTitleOffset(0.30)
 pulls.GetYaxis().SetLabelSize(0.15)
 pulls.GetXaxis().SetTitleSize(0.17)
 pulls.GetXaxis().SetTitleOffset(0.91)
 pulls.GetXaxis().SetLabelSize(0.12)
 pulls.GetXaxis().SetNdivisions(906)
 pulls.GetYaxis().SetNdivisions(305)
 pulls.Draw("same")
 line = ROOT.TLine(1126,0,frame.GetXaxis().GetXmax(),0)
 line1  = ROOT.TLine(1126,1,frame.GetXaxis().GetXmax(),1)
 line2  = ROOT.TLine(1126,-1,frame.GetXaxis().GetXmax(),-1)
 line1.SetLineStyle(2)
 line1.SetLineWidth(2)
 line2.SetLineStyle(2)
 line2.SetLineWidth(2)
 line.Draw("same")
 line1.Draw("same")
 line2.Draw("same")    
 c1.Update()

 canvname+='.pdf'
 c1.SaveAs(canvname)
 c1.SaveAs(canvname.replace("pdf","C"),"C")

def calculateChi2(hdata,nPars,pulls): #THIS NEEDS TO BE FIXED
  
 NumberOfVarBins = 0
 NumberOfObservations_VarBin = 0
 chi2_VarBin = 0.
  
 for p in range (1,pulls.GetNbinsX()+1):
 
     NumberOfVarBins += 1
     #x = ROOT.Double(0.)
     #y = ROOT.Double(0.)
     #pulls.GetPoint(p,x,y)
     #print p,x,y
  
     bin = pulls.GetXaxis().FindBin(x)
     data = pulls.GetBinContent(bin)
     fit_residual = y
     #print "CHI2:",bin,hdata.GetBinCenter(bin),data,fit_residual
     
     if (data>0):
       NumberOfObservations_VarBin+=1
       chi2_VarBin += pow(fit_residual,2)
       
 ndf_VarBin = NumberOfObservations_VarBin - nPars -1 #ndof   
 return [chi2_VarBin,ndf_VarBin] 
 
def makeData(options, dataFile, q, iq, quantiles, hdata, minMJJ=0, maxMJJ=1e+04):
 
 file = h5py.File(options.inputDir+"/"+dataFile,'r')
 # import ipdb; ipdb.set_trace()
 sel_key_q = 'sel_q90' if q == 'q100' else 'sel_' + q # selection column for quantile q (use rejected events of q90 for q100)
 print "Current quantile file: %s, reading quantile %s" % (file, sel_key_q)

 mjj_idx = np.where(file['eventFeatureNames'][()] == 'mJJ')[0]
 sel_idx = np.where(file['eventFeatureNames'][()] == sel_key_q)[0] # 0=rejected 1=accepted
 data = file['eventFeatures'][()] 
 
 if q=='q01':
  for e in range(data.shape[0]):
   #if data[e][mjj_idx] < minMJJ or data[e][mjj_idx] > maxMJJ: continue
   if data[e][sel_idx]==1: hdata.Fill(data[e][mjj_idx])
 elif q=='q100': #if 90% quantile is rejected then events are in the 100-90% slice
  for e in range(data.shape[0]): 
   #if data[e][mjj_idx] < minMJJ or data[e][mjj_idx] > maxMJJ: continue
   if data[e][sel_idx]==0: hdata.Fill(data[e][mjj_idx]) 
 elif q=='total':  
  for e in range(data.shape[0]): hdata.Fill(data[e][mjj_idx]) 
 else:   
  print ".... checking orthogonality wrt",quantiles[iq-1],"quantile...."
  sel_key_iq = 'sel_' + quantiles[iq-1] # selection column for quantile q
  sel_idx_iq = np.where(file['eventFeatureNames'][()] == sel_key_iq)[0] # 0=rejected 1=accepted
  for e in range(data.shape[0]): 
   #if data[e][mjj_idx] < minMJJ or data[e][mjj_idx] > maxMJJ: continue
   if data[e][sel_idx_iq]==0 and data[e][sel_idx]==1: hdata.Fill(data[e][mjj_idx])

def checkSBFit(filename,quantile,roobins,plotname):
 
 fin = ROOT.TFile.Open(filename,'READ')
 workspace = fin.w
 
 model = workspace.pdf('model_s')
 model.Print("v")
 var = workspace.var('mjj')
 data = workspace.data('data_obs')
 
 fres = model.fitTo(data,ROOT.RooFit.SumW2Error(0),ROOT.RooFit.Minos(0),ROOT.RooFit.Verbose(0),ROOT.RooFit.Save(1),ROOT.RooFit.NumCPU(8)) 
 fres.Print()
 
 frame = var.frame()
 data.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data_obs"),ROOT.RooFit.Invisible())
 model.getPdf('JJ_%s'%quantile).plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()), ROOT.RooFit.Binning(roobins))
 model.getPdf('JJ_%s'%quantile).plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model_s"))

 frame3 = var.frame()
 hpull = frame.pullHist("data_obs","model_s",True)
 hpull2 = ROOT.TH1F("hpull2","hpull2",len(binsx)-1, binsx[0], binsx[-1])
 for p in range(hpull.GetN()):
  x = ROOT.Double(0.)
  y = ROOT.Double(0.)
  hpull.GetPoint(p,x,y)
  bin = hpull2.GetXaxis().FindBin(x)
  hpull2.SetBinContent(p+1,y)

 frame3.addPlotable(hpull,"X0 P E1")
 chi2 = frame.chiSquare()
 print chi2
 ndof = 1
 
 data.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data_obs"),ROOT.RooFit.XErrorSize(0))
 #chi2,ndof = calculateChi2(histos_qcd[index],nPars[index],hpull)
 PlotFitResults(frame,fres.GetName(),nPars,frame3,"data_obs","model_s",chi2,ndof,'sbFit_'+plotname)

    
if __name__ == "__main__":

 #python dijetfit.py -i inputdir --sig RSGraviton_WW_NARROW_13TeV_PU40_3.5TeV_parts/RSGraviton_WW_NARROW_13TeV_PU40_3.5TeV_reco.h5 --qcd qcd_sqrtshatTeV_13TeV_PU40_ALL_parts/qcd_sqrtshatTeV_13TeV_PU40_ALL_reco.h5
 #python dijetfit.py -i inputdir --sig RSGraviton_WW_NARROW_13TeV_PU40_1.5TeV_parts/RSGraviton_WW_NARROW_13TeV_PU40_1.5TeV_reco.h5 --qcd qcd_sqrtshatTeV_13TeV_PU40_ALL_parts/qcd_sqrtshatTeV_13TeV_PU40_ALL_reco.h5 --xsec 0.0 -M 1500.0

 #some configuration
 parser = optparse.OptionParser()
 parser.add_option("--xsec","--xsec",dest="xsec",type=float,default=0.0006,help="Injected signal cross section (suggested range 0-0.03)")
 parser.add_option("-M","-M",dest="mass",type=float,default=3500.,help="Injected signal mass")
 parser.add_option("-i","--inputDir",dest="inputDir",default='./',help="directory with all quantiles h5 files")
 parser.add_option("--qcd","--qcd",dest="qcdFile",default='qcd.h5',help="QCD h5 file")
 parser.add_option("--sig","--sig",dest="sigFile",default='signal.h5',help="Signal h5 file")
 parser.add_option("-l","--load_data",dest="load_data",action="store_true",help="Load orthogonal data")
 parser.add_option("--res", "--res", dest="sigRes", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
 (options,args) = parser.parse_args()
  
 xsec = options.xsec
 mass = options.mass
 sigRes = options.sigRes 
 binsx = [1126,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500,5663,5877,6099,6328,6564,6808]
 shift = 1200 - binsx[0] # shift to mjj cut
 binsx = [e+shift for e in binsx]
 roobins = ROOT.RooBinning(len(binsx)-1, array('d',binsx), "mjjbins")
 bins_fine = int(binsx[-1]-binsx[0])
 # quantiles = ['q1','q5','q10','q30','q50','q70','q90','q100','total']
 quantiles = ['q01', 'q10', 'q30', 'q50', 'q70', 'q90','q100','total']
 nPars = 2 # DO THESE NEED TO BE DIFFERENT DEPENDING ON QUANTILE???
 bins_sig_fit = array('f',truncate([binsx[0]+ib for ib in range(bins_fine+1)],0.8*mass,1.2*mass))
 large_bins_sig_fit = array('f',truncate(binsx,0.8*mass,1.2*mass))
 roobins_sig_fit = ROOT.RooBinning(len(large_bins_sig_fit)-1, array('d',large_bins_sig_fit), "mjjbins_sig")

 qr_test_share = 0.8
 qcd_xsec = 8.73e6 # [fb]
 qcd_gen_events = qr_test_share*get_generated_events(options.qcdFile) # 80% of SR qcd events => why applying dEta (SR) cut here but not all the other cuts??
 sig_xsec = 1000. # metric [fb] = 1 [pb]???
 sig_gen_events = qr_test_share*get_generated_events(options.sigFile) # 80% of signal events (but this corresponds to an enormous xsec?!)
 lumi = qcd_gen_events/qcd_xsec # ??? qcd SR lumi
  
 ################################### FIRST PREPARE DATA ###################################
 '''
  MAKE HISTOGRAMS:
  for each quantile + 'bottom rejected' 10% + all events:
    fill mjj histogram for 
    - signal: histos_sig
    - background: histos_qcd
    and save in data_mjj_X_qY.root
 '''       
 histos_sig = []
 histos_qcd = []

 if not options.load_data:
 
  #Signal data preparation 
  for iq,q in enumerate(quantiles):
 
   #histos_sig.append( ROOT.TH1F("mjj_sig_%s"%q,"mjj_sig_%s"%q,bins_fine,binsx[0],binsx[-1]) )
   histos_sig.append( ROOT.TH1F("mjj_sig_%s"%q,"mjj_sig_%s"%q,len(bins_sig_fit)-1,bins_sig_fit) )
   print
   makeData(options,options.sigFile,q,iq,quantiles,histos_sig[-1],0.8*mass,1.2*mass) #first fill orthogonal data histos
   print "************ Found",histos_sig[-1].GetEntries(),"signal events for quantile",q
   print
  
  #Background data preparation
  for iq,q in enumerate(quantiles):
 
   histos_qcd.append( ROOT.TH1F("mjj_qcd_%s"%q,"mjj_qcd_%s"%q,bins_fine,binsx[0],binsx[-1]) )
   print
   makeData(options,options.qcdFile,q,iq,quantiles,histos_qcd[-1]) #first fill orthogonal data histos
   print "************ Found",histos_qcd[-1].GetEntries(),"background events for quantile",q
   print
  
  for h in histos_sig: h.SaveAs("data_"+h.GetName()+".root")
  for h in histos_qcd: h.SaveAs("data_"+h.GetName()+".root")
 
 else: #let's make it faster if you have run once already!
 
  #Load signal data
  for q in quantiles:
  
   fname = "data_mjj_sig_%s.root"%q
   q_datafile = ROOT.TFile.Open(fname,'READ')
   histos_sig.append(q_datafile.Get("mjj_sig_%s"%q))
   histos_sig[-1].SetDirectory(ROOT.gROOT)
   q_datafile.Close()

  #Load background data
  for q in quantiles:
     
   fname = "data_mjj_qcd_%s.root"%q
   q_datafile = ROOT.TFile.Open(fname,'READ')
   histos_qcd.append(q_datafile.Get("mjj_qcd_%s"%q))
   histos_qcd[-1].SetDirectory(ROOT.gROOT)
   q_datafile.Close()
    
  for q,h in enumerate(histos_sig): print "************ Found",h.GetEntries(),"signal events for quantile",quantiles[q]
  for q,h in enumerate(histos_qcd): print "************ Found",h.GetEntries(),"background events for quantile",quantiles[q]

 print "TOTAL SIGNAL EVENTS",histos_sig[-1].GetEntries()
 print "TOTAL BACKGROUND EVENTS",histos_qcd[-1].GetEntries()
 print

 ################################### NOW MAKE THE FITS ###################################
 '''
  for each quantile:
    - fit signal shape -> gauss mu & std ??
    - fit background shape -> exponential lambda ??
    - chi2 ??
 '''

 cmdCombine = 'combineCards.py '
 for iq,q in enumerate(quantiles):
  
  print "########## FIT SIGNAL AND SAVE PARAMETERS ############"
  sig_outfile = ROOT.TFile("sig_fit_%s.root"%q,"RECREATE")
 
  ### create signal model: gaussian centered at mass-center with sigma in {2%,10%of mass center} + crystal ball for asymmetric tail

  fitter=Fitter(['mjj_fine'])
  # if narrow signal, fit narrow model
  if sigRes == "na":
    fitter.signalResonance('model_s',"mjj_fine",mass)
  # else fit broad signal
  else:
    # set crystal ball function params
    alpha_ini=0.85 
    sign_ini=0.5 
    sign_n_stop=20
    # make gauss sigma broader
    sigma = (mass*0.1, mass*0.02, mass*0.2)
    fitter.signalResonance('model_s',"mjj_fine", mass=mass, alpha_ini=alpha_ini, sign_ini=sign_ini, sig_n_stop=sig_n_stop, sigma=sigma)

  ### fit the signal model to sig histogram data

  fitter.w.var("MH").setVal(mass)
  fitter.importBinnedData(histos_sig[iq],['mjj_fine'],'data')
  fres = fitter.fit('model_s','data',[ROOT.RooFit.Save(1)])
  fres.Print()

  ### compute chi-square of compatibility of signal-histogram and signal model for sanity check
  # plot fit result to signal_fit_q.png for each quantile q
 
  mjj_fine = fitter.getVar('mjj_fine')
  mjj_fine.setBins(len(bins_sig_fit))
  chi2_fine = fitter.projection("model_s","data","mjj_fine","signal_fit_%s.png"%q)
  fitter.projection("model_s","data","mjj_fine","signal_fit_%s_log.png"%q,0,True)
  chi2 = fitter.projection("model_s","data","mjj_fine","signal_fit_%s_binned.png"%q,roobins_sig_fit)
  fitter.projection("model_s","data","mjj_fine","signal_fit_%s_log_binned.png"%q,roobins_sig_fit,True)

  # write signal histogram and model with params
 
  sig_outfile.cd()
  histos_sig[iq].Write()
 
  graphs={'mean':ROOT.TGraphErrors(),'sigma':ROOT.TGraphErrors(),'alpha':ROOT.TGraphErrors(),'sign':ROOT.TGraphErrors(),'scalesigma':ROOT.TGraphErrors(),'sigfrac':ROOT.TGraphErrors()}
  for var,graph in graphs.iteritems():
     value,error=fitter.fetch(var)
     graph.SetPoint(0,mass,value)
     graph.SetPointError(0,0.0,error)

  sig_outfile.cd()
  for name,graph in graphs.iteritems(): graph.Write(name)
  
  sig_outfile.Close() 
 
  print "#############################"
  print "signal fit chi2 (fine binning)",chi2_fine
  print "signal fit chi2 (large binning)",chi2
  print "#############################"

  print
  print
  print "############# FIT BACKGROUND AND SAVE PARAMETERS ###########"
  qcd_outfile = ROOT.TFile('qcd_fit_%s.root'%q,'RECREATE')

  ### create signal model: 2-parameter (p1 & p2) exponential
 
  fitter_QCD=Fitter(['mjj_fine'])
  fitter_QCD.qcdShape('model_b','mjj_fine',nPars)
  
  ### fit background model to qcd histogram data

  fitter_QCD.importBinnedData(histos_qcd[iq],['mjj_fine'],'data_qcd') # ??? Fit to qcd histograms (<- we know here that this is qcd?!)
  fres = fitter_QCD.fit('model_b','data_qcd',[ROOT.RooFit.Save(1)])
  fres.Print()
 
  ### compute chi-square of compatibility of qcd-histogram and background model for sanity check
  # plot fit result to qcd_fit_q_binned.png for each quantile q

  chi2_fine = fitter_QCD.projection("model_b","data_qcd","mjj_fine","qcd_fit_%s.png"%q,0,True) # ??? chi2 => sanity check
  chi2_binned = fitter_QCD.projection("model_b","data_qcd","mjj_fine","qcd_fit_%s_binned.png"%q,roobins,True)

  ### write background histogram
 
  qcd_outfile.cd()
  histos_qcd[iq].Write() # ??? => write histo into qcd_outfile

  ### plot background data, fit and 
 
  mjj = fitter_QCD.getVar('mjj_fine')
  mjj.setBins(bins_fine)
  model = fitter_QCD.getFunc('model_b')
  dataset = fitter_QCD.getData('data_qcd')
 
  frame = mjj.frame()
  dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Name("data_qcd"),ROOT.RooFit.Invisible(),ROOT.RooFit.Binning(roobins))
  model.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()),ROOT.RooFit.Binning(roobins))
  model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model_b"),ROOT.RooFit.Binning(roobins))

  framePulls = mjj.frame()
  hpull = frame.pullHist("data_qcd","model_b",True) # ??? pull => second canvas in pull
  #hpull2 = ROOT.TH1F("hpull2","hpull2",len(binsx)-1, binsx[0], binsx[-1])
  #for p in range(hpull.GetN()):
  # x = ROOT.Double(0.)
  # y = ROOT.Double(0.)
  # hpull.GetPoint(p,x,y)
  # #print p,x,y
  # bin = hpull2.GetXaxis().FindBin(x)
  # hpull2.SetBinContent(p+1,y)
 
  framePulls.addPlotable(hpull,"X0 P E1")
  chi2 = frame.chiSquare()
  ndof = 1
  print "chi2 frame:",frame.chiSquare()
 
  dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson),ROOT.RooFit.Name("data_qcd"),ROOT.RooFit.XErrorSize(0),ROOT.RooFit.Binning(roobins))
  #my_chi2,my_ndof = calculateChi2(histos_qcd[index],nPars[index],hpull2)
  #print "my chi2:",chi2/ndof

  # plot full qcd fit results with pull factors to mjj_qcd_q.pdf for each quantile q

  PlotFitResults(frame,fres.GetName(),nPars,framePulls,"data_qcd","model_b",chi2,ndof,histos_qcd[iq].GetName())
 
  # write qcd model with params

  graphs = {}
  for p in range(nPars): graphs['p%i'%(p+1)] = ROOT.TGraphErrors()
  for var,graph in graphs.iteritems():
     print var
     value,error=fitter_QCD.fetch(var)
     graph.SetPoint(0,mass,value)
     graph.SetPointError(0,0.0,error)

  qcd_outfile.cd()          
  for name,graph in graphs.iteritems(): graph.Write(name) # ??? => saving params of bg fit -> load later for datacard
 
  qcd_outfile.Close()

  print "#############################"
  print "bkg fit chi2 (fine binning)",chi2_fine
  print "bkg fit chi2 (large binning)",chi2_binned
  print "bkg fit chi2",chi2
  print "#############################"


  print
  print 
  print "############# GENERATE SIGNAL+BACKGROUND DATA FROM PDFs ###########"
  
  f = ROOT.TFile("/tmp/%s/cache%i.root"%(commands.getoutput("whoami"),random.randint(0, 1e+6)),"RECREATE")
  f.cd()
  w=ROOT.RooWorkspace("w","w")
 
  model_b = fitter_QCD.getFunc('model_b')
  model_s = fitter.getFunc('model_s')
  
  model_b.Print("v")
  model_s.Print("v")
    
  print "Generate",histos_qcd[iq].Integral(),"background events from model_b" # ??? taking qcd number events as they are
  dataqcd = model_b.generateBinned(ROOT.RooArgSet(mjj),histos_qcd[iq].Integral())
  hdataqcd = dataqcd.createHistogram("mjj_fine")
  hdataqcd.SetName("mjj_generate_qcd_%s"%q)

  # signal xsec set to 0 by default, so hdatasig filled with qcd histogram values !

  if xsec != 0:
   print "Generate",int(histos_sig[iq].Integral()*xsec),"signal events from model_s"  # ??? taking sig number events scaled by xsec
   datasig = model_s.generateBinned(ROOT.RooArgSet(mjj),int(histos_sig[iq].Integral()*xsec))
   hdatasig = datasig.createHistogram("mjj_fine")
  else:
   hdatasig = ROOT.TH1F("mjj_generate_sig_%s"%q,"mjj_generate_sig_%s"%q,histos_qcd[iq].GetNbinsX(),histos_qcd[iq].GetXaxis().GetXmin(),histos_qcd[iq].GetXaxis().GetXmax())

  hdatasig.SetName("mjj_generate_sig_%s"%q)
  
  # ??? signal+background fit (total histo)

  sb_outfile = ROOT.TFile('sb_fit_%s.root'%q,'RECREATE')
  sb_outfile.cd()
  htot = ROOT.TH1F()
  htot = hdataqcd.Clone("mjj_generate_tot_%s"%q)
  htot.Add(hdatasig)
  hdatasig.Write("mjj_generate_sig_%s"%q)
  hdataqcd.Write("mjj_generate_qcd_%s"%q)
  htot.Write("mjj_generate_tot_%s"%q)

  w.Delete()
  f.Close()
  f.Delete()
  sb_outfile.Close()
  fitter_QCD.delete()
  fitter.delete()
 
  print
  print
  print "############ MAKE PER CATEGORY DATACARD AND WORKSPACE AND RUN COMBINE #############"
 
  card=DataCardMaker(q)
 
  card.addSignalShape('model_signal_mjj','mjj','sig_fit_%s.root'%q,{'CMS_scale_j':1.0},{'CMS_res_j':1.0})
  # TODO: check numbers!
  # sig_xsec = 1000fb, lumi = qcd SR lumi (although signal in SR AND SB???), sig_gen_events ~ 80% of 900K signal
  # c = number of signal events expected with sig_xsec 1000 fb / number of signals generated = < 1
  constant = sig_xsec*lumi/sig_gen_events  
  if xsec==0: constant = 1.0*constant
  else: constant=xsec*constant
  card.addFixedYieldFromFile('model_signal_mjj',0,'sig_fit_%s.root'%q,histos_sig[iq].GetName(),constant=constant)
  card.addSystematic("CMS_scale_j","param",[0.0,0.012])
  card.addSystematic("CMS_res_j","param",[0.0,0.08]) 
 
  card.addQCDShapeNoTag('model_qcd_mjj','mjj','qcd_fit_%s.root'%q,nPars)
  card.addFloatingYield('model_qcd_mjj',1,'qcd_fit_%s.root'%q,histos_qcd[iq].GetName())
  for i in range(1,nPars+1): card.addSystematic("CMS_JJ_p%i"%i,"flatParam",[])
  card.addSystematic("model_qcd_mjj_JJ_%s_norm"%q,"flatParam",[]) # integral -> anzahl events -> fuer skalierung der genormten roofit histogramm
 
  card.importBinnedData('sb_fit_%s.root'%q,'mjj_generate_tot_%s'%q,["mjj"],'data_obs',1.0)
  card.makeCard()
  card.delete()

  # run combine on datacard -> create workspaces workspace_JJ_0.0_quantile.root
  cmd = 'text2workspace.py datacard_JJ_{label}.txt -o workspace_JJ_{xsec}_{label}.root && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{xsec}_{label} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} --pvalue -n pvalue_{xsec}_{label}'.format(mass=mass,xsec=xsec,label=q)
  print cmd
  os.system(cmd)
  
  #run and visualize s+b fit as sanity check (sb_fit_mjj_qcd_q.root.pdf)
  checkSBFit('workspace_JJ_{xsec}_{label}.root'.format(xsec=xsec,label=q),q,roobins,histos_qcd[iq].GetName()+".root")

 print
 print
 print "############ MAKE N-CATEGORY DATACARD AND WORKSPACE AND RUN COMBINE #############"
 #The difference here is that the background shape comes from one specific quantile (rather than from its own as above)
 for iq,q in enumerate(quantiles):  
 
  if q == 'total': continue
 
  card=DataCardMaker(q+"_4combo")
 
  print "********************** Add signal shape to datacard **********************"
  card.addSignalShape('model_signal_mjj','mjj','sig_fit_%s.root'%q,{'CMS_scale_j':1.0},{'CMS_res_j':1.0})
  # TODO: check data!
  constant = sig_xsec*lumi/sig_gen_events
  if xsec==0: constant = 1.0*constant
  else: constant=xsec*constant
  card.addFixedYieldFromFile('model_signal_mjj',0,'sig_fit_%s.root'%q,histos_sig[iq].GetName(),constant=constant)
  card.addSystematic("CMS_scale_j","param",[0.0,0.012])
  card.addSystematic("CMS_res_j","param",[0.0,0.08]) 
 
  #TAKE BACKGROUND SHAPE COMES FROM BACKGROUND-ENRICHED QUANTILE SLICE --> WHICH ONE? TRY THE Q100 SLICE!
  card.addQCDShapeNoTag('model_qcd_mjj','mjj','qcd_fit_q100.root',nPars) 
  card.addFloatingYield('model_qcd_mjj',1,'qcd_fit_%s.root'%q,histos_qcd[iq].GetName())
  for i in range(1,nPars+1): card.addSystematic("CMS_JJ_p%i"%i,"flatParam",[])
  card.addSystematic("model_qcd_mjj_JJ_q100_4combo_norm","flatParam",[])
 
  card.importBinnedData('sb_fit_%s.root'%q,'mjj_generate_tot_%s'%q,["mjj"],'data_obs',1.0)
  card.makeCard()
  card.delete()
  
  cmdCombine+="quantile_{quantile}=datacard_JJ_{label}.txt ".format(quantile=q,label=q+"_4combo",xsec=xsec)
 
 #MAKE FINAL DATACARD (needs some cosmetics as below) 
 cmdCombine+= '&> datacard_{xsec}_final.txt'.format(xsec=xsec)
 print cmdCombine
 os.system(cmdCombine)
 d = open('datacard_tmp.txt','w')
 dorig = open('datacard_{xsec}_final.txt'.format(xsec=xsec),'r')
 for l in dorig.readlines(): d.write(l)
 d.write('quantile_q100_rate     rateParam       quantile_q100  model_qcd_mjj   1\n')
 d.write('quantile_q90_rate      rateParam       quantile_q90  model_qcd_mjj   (0.20*@0)/0.10  quantile_q100_rate\n')
 d.write('quantile_q70_rate      rateParam       quantile_q70  model_qcd_mjj   (0.20*@0)/0.10  quantile_q100_rate\n')
 d.write('quantile_q50_rate      rateParam       quantile_q50  model_qcd_mjj   (0.20*@0)/0.10  quantile_q100_rate\n')
 d.write('quantile_q30_rate      rateParam       quantile_q30  model_qcd_mjj   (0.20*@0)/0.10  quantile_q100_rate\n') 
 d.write('quantile_q10_rate      rateParam       quantile_q10  model_qcd_mjj   (0.05*@0)/0.10  quantile_q100_rate\n')
 # d.write('quantile_q5_rate      rateParam       quantile_q5  model_qcd_mjj   (0.04*@0)/0.10  quantile_q100_rate\n')
 d.write('quantile_q01_rate      rateParam       quantile_q01  model_qcd_mjj   (0.01*@0)/0.10  quantile_q100_rate\n')
 d.close()
 cmd = 'mv datacard_tmp.txt datacard_{xsec}_final.txt && text2workspace.py datacard_{xsec}_final.txt -o workspace_{xsec}_final.root && combine -M Significance workspace_{xsec}_final.root -m {mass} -n significance_{xsec} && combine -M Significance workspace_{xsec}_final.root -m {mass} --pvalue -n pvalue_{xsec}'.format(mass=mass,xsec=xsec)
 print cmd
 os.system(cmd)
 
 #checkSBFit('workspace_{xsec}_{label}.root'.format(xsec=xsec,label='final'),q+"_4combo",roobins,'final.root')
