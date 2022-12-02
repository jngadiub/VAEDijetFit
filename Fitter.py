import ROOT
import json
from numpy import random
from array import array
import sys,commands
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

class Fitter(object):
    def __init__(self,poi = ['x']):
        self.cache=ROOT.TFile("/tmp/%s/cache%i.root"%(commands.getoutput("whoami"),random.randint(0, 1e+6)),"RECREATE")
        self.cache.cd()

        self.w=ROOT.RooWorkspace("w","w")
        self.dimensions = len(poi)
        self.poi=poi
        for v in poi:
            self.w.factory(v+"[1,161]")

    def delete(self):
     if self.w:
      self.w.Delete()
      self.cache.Close()
      self.cache.Delete()

    def importBinnedData(self,histogram,poi = ["x"],name = "data"):
        cList = ROOT.RooArgList()
        for i,p in enumerate(poi):
            cList.add(self.w.var(p))
            if i==0:
                axis=histogram.GetXaxis()
            elif i==1:
                axis=histogram.GetYaxis()
            elif i==2:
                axis=histogram.GetZaxis()
            else:
                print 'Asking for more than 3 D . ROOT doesnt support that, use unbinned data instead'
                return
            mini=axis.GetXmin()
            maxi=axis.GetXmax()
            bins=axis.GetNbins()
            binningx =[]
            for i in range(1,bins+2):
                #v = mmin + i * (mmax-mmin)/float(N)
                binningx.append(axis.GetBinLowEdge(i))           
            self.w.var(p).setMax(maxi)
            self.w.var(p).setMin(mini)            
            #print " set binning "+str(binningx)
            self.w.var(p).setBinning(ROOT.RooBinning(len(binningx)-1,array("d",binningx)))
            #a = self.w.var(p).getBinning()
            #for b in range(0,a.numBins()+1):
                #print a.binLow(b)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        getattr(self.w,'import')(dataHist,ROOT.RooFit.Rename(name))

    def fetch(self,var):
        self.w.var(var).Print()
        print "Fetching value " ,self.w.var(var).getVal()  
        print "Fetching error " ,self.w.var(var).getError()
        return (self.w.var(var).getVal(), self.w.var(var).getError())

    def getFunc(self,model = "model"):
        return self.w.pdf(model)

    def getData(self,data = "data"):
        return self.w.data(data)

    def getVar(self,var = "mjj"):
        return self.w.var(var)

    def getW(self):
        return self.w
                    
    def fit(self,model = "model",data="data",options=[]):
        if len(options)==0:
            fitresults = self.w.pdf(model).fitTo(self.w.data(data))
        if len(options)==1:
            fitresults = self.w.pdf(model).fitTo(self.w.data(data),options[0])      
        if len(options)==2:
            fitresults = self.w.pdf(model).fitTo(self.w.data(data),options[0],options[1])
        if len(options)==3:
            fitresults = self.w.pdf(model).fitTo(self.w.data(data),options[0],options[1],options[2])
        if len(options)==4:
            fitresults = self.w.pdf(model).fitTo(self.w.data(data),options[0],options[1],options[2],options[3])
          
        if fitresults:
            fitresults.Print() 
            f = ROOT.TFile.Open('fitresults.root','RECREATE')
            fitresults.SetName("fitresults")
            fitresults.Write()
            f.Close()   
    
        return fitresults 

    def getLegend(self):
        self.legend = ROOT.TLegend(0.7510112,0.7183362,0.8502143,0.919833)
        self.legend.SetTextSize(0.032)
        self.legend.SetLineColor(0)
        self.legend.SetShadowColor(0)
        self.legend.SetLineStyle(1)
        self.legend.SetLineWidth(1)
        self.legend.SetFillColor(0)
        self.legend.SetFillStyle(0)
        self.legend.SetMargin(0.35)
        return self.legend
    
    def projection(self,model = "model",data="data",poi="x",filename="fit.root",binning=0,logy=False,xtitle='x',mass=1000):
        
        self.frame=self.w.var(poi).frame()
    
        print "Printing workspace: "
        self.w.Print()
    
        f = ROOT.TFile.Open("fitresults.root",'READ')
        if f: fr = f.Get('fitresults')
        else:
             fr = 0
             print "No fit result found (fitresults.root), plotting model only"
    
        if binning:
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning),ROOT.RooFit.Invisible())
            if fr: self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()),ROOT.RooFit.Binning(binning))
            self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.Binning(binning),ROOT.RooFit.LineColor(ROOT.kRed+1))  
        else: 
            self.w.data(data).plotOn(self.frame,ROOT.RooFit.Invisible())
            if fr: self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.VisualizeError(fr,1),ROOT.RooFit.FillColor(ROOT.kRed-7),ROOT.RooFit.LineColor(ROOT.kRed-7),ROOT.RooFit.Name(fr.GetName()))
            self.w.pdf(model).plotOn(self.frame,ROOT.RooFit.LineColor(ROOT.kRed+1))

        if binning: self.w.data(data).plotOn(self.frame,ROOT.RooFit.Binning(binning))
        else: self.w.data(data).plotOn(self.frame)
    
        self.legend = self.getLegend()
        self.legend.AddEntry( self.w.pdf(model)," Full PDF","l")
                        
        self.c=ROOT.TCanvas("c","c")
        if logy:
          self.frame.SetMinimum(1)
          self.frame.SetMaximum(1e+7)
          self.c.SetLogy()
        self.c.cd()
        self.frame.Draw()
        self.frame.GetYaxis().SetTitle('events')
        self.frame.GetXaxis().SetTitle(xtitle)
        self.frame.SetTitle('')
        self.c.Draw()

        self.legend.Draw("same")        
        self.c.SaveAs(filename)
        pullDist = self.frame.pullHist()
        return self.frame.chiSquare()
         
    def signalResonance(self, name = 'model',poi="MVV",mass=0, sigma=None, alpha=None, sign=None): # poi= mjj_fine
    
        # set default gauss sigma and crystal ball alpha and n tuple if none passed
        sigma = sigma or (mass*0.05, mass*0.02, mass*0.10)
        alpha = alpha or (0.85, 0.60, 1.20)
        sign = sign or (6., 0.1, 150.)

        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    
        # variable[value], [start,stop], [ini,start,stop]
        self.w.factory("MH[1000]")
        # !!! => adapt sigma for broad signal
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(sigma))
        self.w.factory("alpha[%.2f,%.2f,%.2f]"%(alpha))
        self.w.factory("sign[%.1f,%.1f,%.1f]"%(sign))
        self.w.factory("scalesigma[2.0,1.2,3.6]")
        gsigma = ROOT.RooFormulaVar("gsigma","@0*@1", ROOT.RooArgList(self.w.var("sigma"),self.w.var("scalesigma")))
        getattr(self.w,'import')(gsigma,ROOT.RooFit.Rename('gsigma'))
        self.w.factory("Gaussian::gauss(%s,mean,gsigma)"%poi)
        # !!! => adapt crystal ball for asymmetric shapes to broad signal
        self.w.factory("CBShape::cb(%s,mean,sigma,alpha,sign)"%poi) # crystal ball
        self.w.factory('SUM::'+name+'(sigfrac[0.0,0.0,0.850]*gauss,cb)') # sigfrac???

    def signalResonanceDCB(self, name = 'model',poi="MVV",mass=0, singleSided=False):

        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
        self.w.factory("MH[1000]")
        self.w.factory("mean[%.1f,%.1f,%.1f]"%(mass,0.8*mass,1.2*mass))
        self.w.factory("sigma[%.1f,%.1f,%.1f]"%(mass*0.05,mass*0.02,mass*0.10))
        self.w.factory("alpha1[1.8,0.0,2]")
        #self.w.factory("n1[5,0,600]")
        self.w.factory("n1[7]")
        if singleSided:
            self.w.factory("alpha2[1000000.0]")
            self.w.factory("n2[0]")
        else:
            self.w.factory("alpha2[1.2,0.0,10]")
            #self.w.factory("n2[5,0,50]")
            self.w.factory("n2[4]")
        peak_vv = ROOT.RooDoubleCB(name,'modelS',self.w.var(poi),self.w.var('mean'),self.w.function('sigma'),self.w.var('alpha1'),self.w.var('n1'),self.w.var('alpha2'),self.w.var('n2'))
        getattr(self.w,'import')(peak_vv,ROOT.RooFit.Rename(name))

    def qcdShape(self,name = 'model',poi="MVV",nPars=2):
    
        ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
    
        if nPars==2:
         self.w.factory("p1[9.28433e+00, -100. , 100.]") # ??? where do these values come from? => from cms diboson analysis
         self.w.factory("p2[1.03641e+01, -200, 200]")
         model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2")))   
        elif nPars==3: 
         self.w.factory("p1[9.28433e+00, -100. , 100.]")
         self.w.factory("p2[1.03641e+01, -200, 200]")
         self.w.factory("p3[2.35256e+00, -100., 100.]")
         model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3")))
        elif nPars==4:
         #default set
         self.w.factory("p1[9.28433e+00, -100. , 100.]")
         self.w.factory("p2[1.03641e+01, -200, 200]")    
         self.w.factory("p3[2.35256e+00, -100., 100.]")
         self.w.factory("p4[4.17695e-01, -100., 100.]")
         #set-1
         #self.w.factory("p1[3.64825e-05, -1000. , 1000.]")
         #self.w.factory("p2[2.78348e+00, -100, 100]")   
         #self.w.factory("p3[7.17321e+00, 0., 50.]")
         #self.w.factory("p4[-9.83737e-01, -1000., 1000.]")
         #set-2
         #self.w.factory("p1[3.45939e+02, -100. , 200.]")
         #self.w.factory("p2[-4.45967e+00, -100, 100]")  
         #self.w.factory("p3[1.31413e+00, 0., 50.]")
         #self.w.factory("p4[1.91782e+02, -1000., 1000.]")
         #set-3
         #self.w.factory("p1[4.27857e+01, -100. , 1000.]")
         #self.w.factory("p2[-1.37580e-02, -100, 100]")  
         #self.w.factory("p3[7.46099e+00, 0., 50.]")
         #self.w.factory("p4[1.78503e+00, -1000., 100.]")
         model = ROOT.RooGenericPdf(name, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4")))
         #alt func
         #model = ROOT.RooGenericPdf(name, "( @1*pow(1-@0/13000 + @4*pow(@0/13000,2),@2) ) / ( pow(@0/13000,@3) )", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4")))
        elif nPars ==5:
            self.w.factory("p1[9.28433e+00, -100. , 100.]")
            self.w.factory("p2[1.03641e+01, -200, 200]")	 
            self.w.factory("p3[2.35256e+00, -100., 100.]")
            self.w.factory("p4[4.17695e-01, -100., 100.]")
            self.w.factory("p5[1.00000e+01, -100., 100.]")
            model = ROOT.RooGenericPdf(name, "pow(exp(-@0/13000.),@4) *pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4"), self.w.var("p5")))
        elif nPars==6:
            self.w.factory("p1[6.02992e+00, -100. , 100]")
            self.w.factory("p2[6.28634e+00, -200., 200]")
            self.w.factory("p3[6.32552e-01, -100., 100.]")
            self.w.factory("p4[-8.53977e-02, -100., 100.]")
            self.w.factory("p5[1000, 1., 5000.]")
            self.w.factory("p6[1600, 0., 2500.]")
            model = ROOT.RooGenericPdf(name, "(0.5*tanh((@0-@6)/@5) + .5)*pow(1-@0/13000., @1) / ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(poi), self.w.var("p1"), self.w.var("p2"), self.w.var("p3"), self.w.var("p4"), self.w.var("p5"), self.w.var("p6"))) 

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(name))
