import h5py, math, commands, random
from array import array
import numpy as np
import time, sys, os, optparse

import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)

from Fitter import Fitter
from DataCardMaker import DataCardMaker
        			        		       
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
 
 
if __name__ == "__main__":

 #some configuration
 parser = optparse.OptionParser()
 parser.add_option("--xsec","--xsec",dest="xsec",type=float,default=0.025,help="Injected signal cross section (suggested range 0-0.03)")
 parser.add_option("-M","-M",dest="mass",type=float,default=3500.,help="Injected signal mass")
 parser.add_option("--index","--index",dest="index",type=int,default=1,help="Which selections to fit (0=inclusive ; 1=accepted ; 2=rejected)")
 parser.add_option("--qcd","--qcd",dest="qcdFile",default='qcd.h5',help="QCD h5 file")
 parser.add_option("--sig","--sig",dest="sigFile",default='signal.h5',help="Signal h5 file")
 parser.add_option("--twoCatFit","--twoCatFit",dest="twoCatFit",action="store_true",help="Run two categories fit")
 (options,args) = parser.parse_args()
 
 xsec = options.xsec
 index = options.index 
 mass = options.mass 
 labels=['total','accepted','rejected']
 binsx = [1126,1181,1246,1313,1383,1455,1530,1607,1687,1770,1856,1945,2037,2132,2231,2332,2438,2546,2659,2775,2895,3019,3147,3279,3416,3558,3704,3854,4010,4171,4337,4509,4686,4869,5058,5253,5500,5663,5877,6099,6328,6564,6808]
 roobins = ROOT.RooBinning(len(binsx)-1, array('d',binsx), "mjjbins")
 bins_fine = int(binsx[-1]-binsx[0])
 
 
 #make signal histograms with fine binning
 h5FileSig = h5py.File(options.sigFile,'r') 

 histos_sig = []
 for l in labels: histos_sig.append( ROOT.TH1F("mjj_sig_%s"%l,"mjj_sig_%s"%l,bins_fine,binsx[0],binsx[-1]) )
 for h in histos_sig: h.SetBinErrorOption(ROOT.TH1.kPoisson)
 
 mjj_idx = np.where(h5FileSig['eventFeatureNames'][()] == 'mJJ')[0]
 sel_idx = np.where(h5FileSig['eventFeatureNames'][()] == 'sel')[0] # 0=rejected 1=accepted
 sig = h5FileSig['eventFeatures'][()]
 for e in range(sig.shape[0]):
  if sig[e][sel_idx] == 0: histos_sig[2].Fill(sig[e][mjj_idx]) #rejected
  if sig[e][sel_idx] == 1: histos_sig[1].Fill(sig[e][mjj_idx]) #accepted
 
 histos_sig[0].Add(histos_sig[1])
 histos_sig[0].Add(histos_sig[2])   
 
 h5FileSig.close()

 
 #make qcd histograms with fine binning
 h5FileBkg = h5py.File(options.qcdFile,'r') 

 histos_qcd = []
 for l in labels:
  histos_qcd.append( ROOT.TH1F("mjj_qcd_%s"%l,"mjj_qcd_%s"%l,bins_fine,binsx[0],binsx[-1]) )
  #histos_qcd.append(ROOT.TH1F("mjj_qcd_%s"%l,"mjj_qcd_%s"%l,len(binsx)-1,array('f',binsx)))
 for h in histos_qcd: h.SetBinErrorOption(ROOT.TH1.kPoisson)
 
 mjj_idx = np.where(h5FileBkg['eventFeatureNames'][()] == 'mJJ')[0]
 sel_idx = np.where(h5FileBkg['eventFeatureNames'][()] == 'sel')[0] # 0=rejected 1=accepted
 bkg = h5FileBkg['eventFeatures'][()]
 for e in range(bkg.shape[0]):
  if bkg[e][sel_idx] == 0: histos_qcd[2].Fill(bkg[e][mjj_idx]) #rejected
  if bkg[e][sel_idx] == 1: histos_qcd[1].Fill(bkg[e][mjj_idx]) #accepted
 
 histos_qcd[0].Add(histos_qcd[1])
 histos_qcd[0].Add(histos_qcd[2])   

 h5FileBkg.close()
      
 ########## FIT SIGNAL AND SAVE PARAMETERS ############
 sig_outfile = ROOT.TFile("sig_fit_%s.root"%labels[index],"RECREATE")
 
 fitter=Fitter(['mjj_fine'])
 fitter.signalResonance('model_s',"mjj_fine",mass)
 fitter.w.var("MH").setVal(mass)
 fitter.importBinnedData(histos_sig[index],['mjj_fine'],'data')
 fres = fitter.fit('model_s','data',[ROOT.RooFit.Save(1)])
 fres.Print()
 
 mjj_fine = fitter.getVar('mjj_fine')
 mjj_fine.setBins(bins_fine)
 chi2_fine = fitter.projection("model_s","data","mjj_fine","signal_fit_%s.png"%labels[index])

 fitter.projection("model_s","data","mjj_fine","signal_fit_%s_log.png"%labels[index],0,True)
 chi2 = fitter.projection("model_s","data","mjj_fine","signal_fit_%s_binned.png"%labels[index],roobins)
 fitter.projection("model_s","data","mjj_fine","signal_fit_%s_log_binned.png"%labels[index],roobins,True)
 
 sig_outfile.cd()
 histos_sig[index].Write()
 
 graphs={'mean':ROOT.TGraphErrors(),'sigma':ROOT.TGraphErrors(),'alpha':ROOT.TGraphErrors(),'sign':ROOT.TGraphErrors(),'scalesigma':ROOT.TGraphErrors(),'sigfrac':ROOT.TGraphErrors()}
 for var,graph in graphs.iteritems():
     value,error=fitter.fetch(var)
     graph.SetPoint(0,mass,value)
     graph.SetPointError(0,0.0,error)

 sig_outfile.cd()
 for name,graph in graphs.iteritems():
    graph.Write(name)
  
 sig_outfile.Close() 
 
 print "#############################"
 print "signal fit chi2 (fine binning)",chi2_fine
 print "signal fit chi2 (large binning)",chi2
 print "#############################"
 time.sleep(20) 
  
 ############# FIT QCD AND SAVE PARAMETERS ###########
 qcd_outfile = ROOT.TFile('qcd_fit_%s.root'%labels[index],'RECREATE')
  
 nPars = [2,2,2]
 
 fitter_QCD=Fitter(['mjj_fine'])
 fitter_QCD.qcdShape('model_b','mjj_fine',nPars[index])
 fitter_QCD.importBinnedData(histos_qcd[index],['mjj_fine'],'data_qcd')
 fres = fitter_QCD.fit('model_b','data_qcd',[ROOT.RooFit.Save(1)])
 fres.Print()
 
 chi2_fine = fitter_QCD.projection("model_b","data_qcd","mjj_fine","qcd_fit_%s.png"%labels[index],0,True)
 chi2_binned = fitter_QCD.projection("model_b","data_qcd","mjj_fine","qcd_fit_%s_binned.png"%labels[index],roobins,True)
 
 qcd_outfile.cd()
 histos_qcd[index].Write()
 
 mjj = fitter_QCD.getVar('mjj_fine')
 mjj.setBins(bins_fine)
 model = fitter_QCD.getFunc('model_b')
 dataset = fitter_QCD.getData('data_qcd')
 
 frame = mjj.frame()
 dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Name("data_qcd"),ROOT.RooFit.Invisible(),ROOT.RooFit.Binning(roobins))
 model.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()),ROOT.RooFit.Binning(roobins))
 model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model_b"),ROOT.RooFit.Binning(roobins))

 framePulls = mjj.frame()
 hpull = frame.pullHist("data_qcd","model_b",True)
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
 PlotFitResults(frame,fres.GetName(),nPars[index],framePulls,"data_qcd","model_b",chi2,ndof,histos_qcd[index].GetName())
 time.sleep(10)
 
 graphs = {}
 for p in range(nPars[index]): graphs['p%i'%(p+1)] = ROOT.TGraphErrors()
 for var,graph in graphs.iteritems():
     print var
     value,error=fitter_QCD.fetch(var)
     graph.SetPoint(0,mass,value)
     graph.SetPointError(0,0.0,error)

 qcd_outfile.cd()          
 for name,graph in graphs.iteritems():
    graph.Write(name)
 
 qcd_outfile.Close()

 print "#############################"
 print "bkg fit chi2 (fine binning)",chi2_fine
 print "bkg fit chi2 (large binning)",chi2_binned
 print "bkg fit chi2",chi2
 print "#############################"
 time.sleep(20) 
   
 ############# SIGNAL+BACKGROUND FIT ###########
 fin_qcd = ROOT.TFile.Open('qcd_fit_%s.root'%labels[index],'READ')
 qcd_data_th1 = fin_qcd.Get(histos_qcd[index].GetName()) 
 
 sb_outfile = ROOT.TFile('sb_fit_%s.root'%labels[index],'RECREATE')
 sb_outfile.cd()
 htot = ROOT.TH1F()
 htot = qcd_data_th1.Clone("htot")
 hsig_generate = ROOT.TH1F("mjj_generate_sig_%s"%l,"mjj_generate_sig_%s"%l,bins_fine,binsx[0],binsx[-1])
 hsig_generate.FillRandom(histos_sig[index],int(histos_sig[index].Integral()*xsec))
 htot.Add(hsig_generate)
 htot.Write()
 qcd_data_th1.Write('qcd')
 hsig_generate.Write('signal')
    
 f = ROOT.TFile("/tmp/%s/cache%i.root"%(commands.getoutput("whoami"),random.randint(0, 1e+6)),"RECREATE")
 f.cd()
 w=ROOT.RooWorkspace("w","w")
 
 model_b = fitter_QCD.getFunc('model_b')
 model_s = fitter.getFunc('model_s')
  
 model_b.Print("v")
 model_s.Print("v")

 Ns = ROOT.RooRealVar("Ns", "signal yield",hsig_generate.Integral(),hsig_generate.Integral()-2*math.sqrt(hsig_generate.Integral()),hsig_generate.Integral()+2*math.sqrt(hsig_generate.Integral()))
 Nb = ROOT.RooRealVar("Nb", "background yield", histos_qcd[index].Integral(), 0, histos_qcd[index].Integral()*2.)
 
 model = ROOT.RooAddPdf("model", "gaussian plus exponential PDF", ROOT.RooArgList(model_s, model_b), ROOT.RooArgList(Ns, Nb))
 getattr(w,'import')(model)

 dataset = ROOT.RooDataHist("data_obs", "data_obs", ROOT.RooArgList(mjj_fine), ROOT.RooFit.Import(htot))  
 getattr(w,'import')(dataset)

 fres = model.fitTo(dataset,ROOT.RooFit.Save())
 fres.Print()

 frame = mjj_fine.frame()
 dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data_obs"),ROOT.RooFit.Invisible())
 model.plotOn(frame,ROOT.RooFit.VisualizeError(fres,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fres.GetName()), ROOT.RooFit.Binning(roobins))
 model.plotOn(frame,ROOT.RooFit.LineColor(ROOT.kRed+1),ROOT.RooFit.Name("model"))

 frame3 = mjj.frame()
 hpull = frame.pullHist("data_obs","model",True)
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
 
 dataset.plotOn(frame,ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), ROOT.RooFit.Binning(roobins),ROOT.RooFit.Name("data"),ROOT.RooFit.XErrorSize(0))
 #chi2,ndof = calculateChi2(histos_qcd[index],nPars[index],hpull)
 PlotFitResults(frame,fres.GetName(),nPars[index],frame3,"data_obs","model",chi2,ndof,'sbFit_'+histos_qcd[index].GetName())

 w.Delete()
 f.Close()
 f.Delete()
 sb_outfile.Close()
 fitter_QCD.delete()
 fitter.delete()
   
 ############ MAKE DATACARD AND WORKSPACE AND RUN COMBINE #############
 
 card=DataCardMaker(labels[index])
 
 if not options.twoCatFit or (options.twoCatFit and index!=2):
  print "********************** Add signal shape to datacard **********************"
  card.addSignalShape('model_signal_mjj','mjj','sig_fit_%s.root'%labels[index],{'CMS_scale_j':1.0},{'CMS_res_j':1.0})
  if xsec==0: constant = 1.0
  else: constant=xsec
  card.addFixedYieldFromFile('model_signal_mjj',0,'sig_fit_%s.root'%labels[index],histos_sig[index].GetName(),constant=constant)
  #card.addSystematic("CMS_scale_j","param",[0.0,1.0])
  #card.addSystematic("CMS_res_j","param",[0.0,1.0]) 
 
 card.addQCDShapeNoTag('model_qcd_mjj','mjj','qcd_fit_%s.root'%labels[index],nPars[index])
 card.addFloatingYield('model_qcd_mjj',1,'qcd_fit_%s.root'%labels[index],histos_qcd[index].GetName())
 for i in range(1,nPars[index]+1): card.addSystematic("CMS_JJ_p%i"%i,"flatParam",[])
 if not options.twoCatFit or (options.twoCatFit and index==2): card.addSystematic("model_qcd_mjj_JJ_%s_norm"%labels[index],"flatParam",[])
 
 if options.twoCatFit and index==2: card.importBinnedData('sb_fit_%s.root'%labels[index],'qcd',["mjj"],'data_obs',1.0)
 else: card.importBinnedData('sb_fit_%s.root'%labels[index],'htot',["mjj"],'data_obs',1.0)
 card.makeCard()
 card.delete()

 if not options.twoCatFit:
  cmd = 'text2workspace.py datacard_JJ_{label}.txt -o workspace_JJ_{label}.root && combine -M Significance workspace_JJ_{label}.root -m {mass} -n significance_{xsec}_{label} && combine -M Significance workspace_JJ_{label}.root -m {mass} --pvalue -n pvalue_{xsec}_{label}'.format(mass=mass,xsec=xsec,label=labels[index])
  print cmd
  os.system(cmd)
 elif options.twoCatFit and index!=2:
  cmd = 'combineCards.py signal_region=datacard_JJ_accepted.txt control_region=datacard_JJ_rejected.txt &> datacard.txt'
  print cmd
  os.system(cmd)
  d = open('datacard_tmp.txt','w')
  dorig = open('datacard.txt','r')
  for l in dorig.readlines(): d.write(l)
  d.write('qcd_control_region_rate	rateParam	control_region	model_qcd_mjj	1\n')
  d.write('qcd_signal_region_rate	rateParam	signal_region	model_qcd_mjj	(0.05*@0)/0.95	qcd_control_region_rate\n')
  d.close()
  cmd = 'mv datacard_tmp.txt datacard.txt && text2workspace.py datacard.txt -o workspace.root && combine -M Significance workspace.root -m {mass} -n significance_{xsec} && combine -M Significance workspace.root -m {mass} --pvalue -n pvalue_{xsec}'.format(mass=mass,xsec=xsec)
  os.system(cmd)
