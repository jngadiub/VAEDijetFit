import ROOT
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
import json
import sys
import os
from array import array
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

class DataCardMaker:
    def __init__(self,tag, outDir):
        self.systematics=[]
        self.tag="JJ"+"_"+tag
        self.outDir = outDir
        self.rootFile = ROOT.TFile(os.path.join(outDir, "datacardInputs_%s.root"%self.tag),"RECREATE")
        self.rootFile.cd()
        self.w=ROOT.RooWorkspace("w","w")
        self.luminosity = 1.0
        self.contributions=[]
        self.systematics=[]

    def delete(self):
     if self.w:
      self.w.Delete()
      self.rootFile.Close()
      self.rootFile.Delete()
      
    def makeCard(self):

        f = open(os.path.join(self.outDir,"datacard_"+self.tag+'.txt'),'w')
        datacard_inputs_file = "datacardInputs_"+self.tag
        f.write('imax 1\n')
        f.write('jmax {n}\n'.format(n=len(self.contributions)-1))
        f.write('kmax *\n')
        f.write('-------------------------\n')
        for c in self.contributions:
            f.write('shapes {name} {channel} {file}.root w:{pdf}\n'.format(name=c['name'], channel=self.tag, file=datacard_inputs_file, pdf=c['pdf']))
        f.write('shapes {name} {channel} {file}.root w:{name}\n'.format(name="data_obs", channel=self.tag, file=datacard_inputs_file))
        f.write('-------------------------\n')
        f.write('bin '+self.tag+'\n')
        f.write('observation  -1\n')
        f.write('-------------------------\n')
        f.write('bin\t') 

        for shape in self.contributions: f.write(self.tag+'\t')
        f.write('\n')

        #Sort the shapes by ID 
 
        shapes = sorted(self.contributions,key=lambda x: x['ID'])
        #print names
        f.write('process\t')
        for shape in shapes:
            f.write(shape['name']+'\t')
        f.write('\n')

        #Print ID
        f.write('process\t')
        for shape in shapes:
            f.write(str(shape['ID'])+'\t')
        f.write('\n')

        #print rates
        f.write('rate\t')
        for shape in shapes:
            f.write(str(shape['yield'])+'\t')
        f.write('\n')


        #Now systematics
        for syst in self.systematics:
            if syst['kind'] == 'param':
                f.write(syst['name']+'\t'+'param\t' +str(syst['values'][0])+'\t'+str(syst['values'][1])+'\n')
            elif syst['kind'] == 'flatParam':
                f.write(syst['name']+'\t'+'flatParam\n')
            elif 'rateParam' in syst['kind']:
                line = syst['name']+'\t'+str(syst['kind'])+"\t"+str(syst['bin'])+"\t"+str(syst['process'])+"\t"+str(syst['values'])+"\t"+str(syst['variables'])
                line+='\n' 
                f.write(line)
                
            elif syst['kind'] == 'discrete':
                f.write(syst['name']+'\t'+'discrete\n')

            elif syst['kind'] == 'lnN' :
                f.write(syst['name']+'\t'+ 'lnN\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
            elif syst['kind'] == 'lnU': 
                f.write(syst['name']+'\t'+ 'lnU\t' )
                for shape in shapes:
                    has=False
                    for name,v in syst['values'].iteritems():
                        if shape['name']==name:
                            f.write(str(v)+'\t' )
                            has=True
                            break;
                    if not has:
                            f.write('-\t' )
                f.write('\n' )
                            
                        
        f.close()


        self.rootFile.cd()
        self.w.Write("w",0,2000000000)
        self.rootFile.Close()
        
    def importBinnedData(self,filename,histoname,poi,name = "data_obs",scale=1):
        f=ROOT.TFile(filename)
        histogram=f.Get(histoname)
        histogram.Scale(scale)
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
            self.w.var(p).setMin(mini)
            self.w.var(p).setMax(maxi)
            self.w.var(p).setBins(bins)
            self.w.var(p).setBins(bins,"cache")
            mjj=self.w.var(p)
        dataHist=ROOT.RooDataHist(name,name,cList,histogram)
        #dataHist=ROOT.RooDataHist(name, name, ROOT.RooArgList(mjj), ROOT.RooFit.Import(histogram)) 
        getattr(self.w,'import')(dataHist,ROOT.RooFit.RenameVariable(name,name))
        
    def addSystematic(self,name,kind,values,bin="",process="",variables="",addPar = ""):
        if kind != 'rateParam': self.systematics.append({'name':name,'kind':kind,'values':values })
        else: self.systematics.append({'name':name,'kind':kind,'bin':bin,'process':process,'values':values,'variables':variables})
        
    def addFixedYieldFromFile(self,name,ID,filename,histoName,constant=1.0):
        pdfName="_".join([name,self.tag])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.GetEntries()*self.luminosity*constant # !!!
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':events})

    # add a floatable number of events value
    def addFloatingYield(self,name,ID,filename,histoName,mini=0,maxi=1e+7,constant=False):
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.GetEntries()
        self.w.factory("{name}[{val},{mini},{maxi}]".format(name=pdfNorm,val=events,mini=mini,maxi=maxi))       
        if constant:
            self.w.var(pdfNorm).setConstant(1)
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})

    def addFloatingYieldCorr(self,name,ID,filename,histoName,fraction,constant=1.0,mini=0,maxi=1e+7):
        pdfName="_".join([name,self.tag])
        pdfNorm="_".join([name,self.tag,"norm"])
        f=ROOT.TFile(filename)
        histogram=f.Get(histoName)
        events=histogram.GetEntries()
        self.w.factory("{name}[{val},{mini},{maxi}]".format(name="shapeBkg_model_qcd_mjj_JJ_q0__norm",val=constant,mini=mini,maxi=maxi)) #the value here can be whatever
        fName = "_".join([name,self.tag,"fraction"])
        self.w.factory("{name}[{val},{mini},{maxi}]".format(name=fName,val=fraction,mini=0,maxi=1))
        self.w.var(fName).setConstant(1)
        prod = ROOT.RooFormulaVar(pdfNorm,"@0*@1", ROOT.RooArgList(self.w.var("shapeBkg_model_qcd_mjj_JJ_q0__norm"),self.w.var(fName)))
        getattr(self.w,'import')(prod,ROOT.RooFit.Rename(pdfNorm))
        self.contributions.append({'name':name,'pdf':pdfName,'ID':ID,'yield':1.0})

    def addSignalShape(self, name, variable, jsonFile, scale ={}, resolution={}):
    
        pdfName="_".join([name,self.tag])
        
        #self.w.factory("MH[3000]")
        #self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
       
        self.w.factory(variable+"[0,13000]")
        
        f = ROOT.TFile(jsonFile,'READ')
        meanG = f.Get('mean')
        sigmaG = f.Get('sigma')
        alphaG = f.Get('alpha')
        scalesigmaG = f.Get('scalesigma')
        sigfracG = f.Get('sigfrac')
        signG = f.Get('sign')
        
        x = ROOT.Double(0.)
        mean = ROOT.Double(0.)
        meanG.GetPoint(0,x,mean)
        sigma = ROOT.Double(0.)
        sigmaG.GetPoint(0,x,sigma)      
        alpha = ROOT.Double(0.)
        alphaG.GetPoint(0,x,alpha)
        scalesigma = ROOT.Double(0.)
        scalesigmaG.GetPoint(0,x,scalesigma)
        sigfrac = ROOT.Double(0.)
        sigfracG.GetPoint(0,x,sigfrac)
        sign = ROOT.Double(0.)
        signG.GetPoint(0,x,sign)
                
        
        meanVar = "_".join(["MEAN",name,self.tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=meanVar,param=mean,vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))

        sigmaVar = "_".join(["SIGMA",name,self.tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=sigmaVar,param=sigma,vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
                
        alphaVar = "_".join(["ALPHA",name,self.tag])            
        alpha = ROOT.RooRealVar(alphaVar,alphaVar,alpha)
        getattr(self.w,'import')(alpha,ROOT.RooFit.Rename(alphaVar))
        
        sigfracVar = "_".join(["SIGFRAC",name,self.tag])
        sigfrac = ROOT.RooRealVar(sigfracVar,sigfracVar,sigfrac)
        getattr(self.w,'import')(sigfrac,ROOT.RooFit.Rename(sigfracVar))
        
        scalesigmaVar = "_".join(["SCALESIGMA",name,self.tag])
        scalesigma = ROOT.RooRealVar(scalesigmaVar,scalesigmaVar,scalesigma)
        getattr(self.w,'import')(scalesigma,ROOT.RooFit.Rename(scalesigmaVar))
        
        signVar = "_".join(["SIGN",name,self.tag])
        sign = ROOT.RooRealVar(signVar,signVar,sign)    
        getattr(self.w,'import')(sign,ROOT.RooFit.Rename(signVar))

        gsigmaVar = "_".join(["GSIGMA",name,self.tag])          
        gsigma = ROOT.RooFormulaVar(gsigmaVar,"@0*@1", ROOT.RooArgList(self.w.function(sigmaVar),scalesigma))
        #getattr(self.w,'import')(gsigma,ROOT.RooFit.Rename(gsigmaVar))      

        gaussFunc = "_".join(["gauss",name,self.tag])   
        gauss = ROOT.RooGaussian(gaussFunc, gaussFunc, self.w.var(variable), self.w.function(meanVar), gsigma)
        cbFunc = "_".join(["cb",name,self.tag])
        cb    = ROOT.RooCBShape(cbFunc, cbFunc,self.w.var(variable), self.w.function(meanVar), self.w.function(sigmaVar), alpha, sign)
        model = ROOT.RooAddPdf(pdfName, pdfName, gauss, cb, self.w.var(sigfracVar))     
        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))

    def addSignalShapeDCB(self, name, variable, jsonFile, scale ={}, resolution={}):
    
        pdfName="_".join([name,self.tag])
        
        #self.w.factory("MH[3000]")
        #self.w.var("MH").setConstant(1)
       
        scaleStr='0'
        resolutionStr='0'

        scaleSysts=[]
        resolutionSysts=[]
        for syst,factor in scale.iteritems():
            self.w.factory(syst+"[0,-0.1,0.1]")
            scaleStr=scaleStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            scaleSysts.append(syst)
        for syst,factor in resolution.iteritems():
            self.w.factory(syst+"[0,-0.5,0.5]")
            resolutionStr=resolutionStr+"+{factor}*{syst}".format(factor=factor,syst=syst)
            resolutionSysts.append(syst)
       
        self.w.factory(variable+"[0,13000]")
        
        f = ROOT.TFile(jsonFile,'READ')
        G_mean = f.Get('mean')
        G_sigma = f.Get('sigma')
        G_alpha1 = f.Get('alpha1')
        G_alpha2 = f.Get('alpha2')
        G_n1 = f.Get('n1')
        G_n2 = f.Get('n2')
        
        x = ROOT.Double(0.)
        mean = ROOT.Double(0.)
        G_mean.GetPoint(0,x,mean)
        sigma = ROOT.Double(0.)
        G_sigma.GetPoint(0,x,sigma)      
        alpha1 = ROOT.Double(0.)
        G_alpha1.GetPoint(0,x,alpha1)
        alpha2 = ROOT.Double(0.)
        G_alpha2.GetPoint(0,x,alpha2)
        n1 = ROOT.Double(0.)
        G_n1.GetPoint(0,x,n1)                
        n2 = ROOT.Double(0.)
        G_n2.GetPoint(0,x,n2)  

        meanVar = "_".join(["MEAN",name,self.tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=meanVar,param=mean,vv_syst=scaleStr,vv_systs=','.join(scaleSysts)))

        sigmaVar = "_".join(["SIGMA",name,self.tag])
        self.w.factory("expr::{name}('{param}*(1+{vv_syst})',{vv_systs},{param})".format(name=sigmaVar,param=sigma,vv_syst=resolutionStr,vv_systs=','.join(resolutionSysts)))
                
        alpha1Var = "_".join(["ALPHA1",name,self.tag])            
        alpha1 = ROOT.RooRealVar(alpha1Var,alpha1Var,alpha1)
        getattr(self.w,'import')(alpha1,ROOT.RooFit.Rename(alpha1Var))
        
        alpha2Var = "_".join(["ALPHA2",name,self.tag])            
        alpha2 = ROOT.RooRealVar(alpha2Var,alpha2Var,alpha2)
        getattr(self.w,'import')(alpha2,ROOT.RooFit.Rename(alpha2Var))
        
        n1Var = "_".join(["N1",name,self.tag])
        n1 = ROOT.RooRealVar(n1Var,n1Var,n1)
        getattr(self.w,'import')(n1,ROOT.RooFit.Rename(n1Var))

        n2Var = "_".join(["N2",name,self.tag])
        n2 = ROOT.RooRealVar(n2Var,n2Var,n2)
        getattr(self.w,'import')(n2,ROOT.RooFit.Rename(n2Var))        
    
        pdfName="_".join([name,self.tag])
        vvMass = ROOT.RooDoubleCB(pdfName,pdfName,self.w.var(variable),self.w.function(meanVar),self.w.function(sigmaVar),self.w.function(alpha1Var),self.w.function(n1Var),self.w.function(alpha2Var),self.w.function(n2Var))
        getattr(self.w,'import')(vvMass,ROOT.RooFit.RenameVariable(pdfName,pdfName))

    def addQCDShape(self,name,variable,preconstrains,nPars=4):

        pdfName=name+"_"+self.tag
        
        MVV=variable
        if self.w.var(MVV) == None: self.w.factory(MVV+"[0,10000]")
        
        errs = [100,200,100,100,100]
        values = [9.28433e+00,1.03641e+01,2.35256e+00,4.17695e-01,1.00000e+01]
        f = ROOT.TFile.Open(preconstrains,'READ')
        parsG = [f.Get('p%i'%i) for i in range(1,nPars+1)]
        pars_val = [ROOT.Double(0.) for i in range(0,nPars)]       
        for i in range(1,nPars+1):
         x = ROOT.Double(0.)         
         pName="_".join(["CMS_JJ_p%i"%i,self.tag])
         if nPars==6 and i==5:
            errUp=5000
            errDown=1
         #   pars_val[i-1] = 1000
         elif nPars==6 and i==6:
            errUp=2500
            errDown=0
         #   pars_val[i-1] = 1600
         else:
          errUp=errs[i-1]
          errDown=-errs[i-1]
         #pars_val[i-1] = values[i-1]
         parsG[i-1].GetPoint(0,x,pars_val[i-1])
         #errUp = pars_val[i-1]+parsG[i-1].GetErrorYhigh(0)*1000.
         #errDown = -pars_val[i-1]-parsG[i-1].GetErrorYlow(0)*1000.
         print i,pName,pars_val[i-1],parsG[i-1].GetErrorYhigh(0),parsG[i-1].GetErrorYlow(0),errUp,errDown
         self.w.factory("{name}[{val},{errDown},{errUp}]".format(name=pName,val=pars_val[i-1],errUp=errUp,errDown=errDown))
        
        if nPars==2: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag)))  
        elif nPars==3: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag)))
        elif nPars==4: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag), self.w.var("CMS_JJ_p4_%s"%self.tag)))
        elif nPars==5: model = ROOT.RooGenericPdf(pdfName, "pow(exp(-@0/13000.),@4) * pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag), self.w.var("CMS_JJ_p4_%s"%self.tag), self.w.var("CMS_JJ_p5_%s"%self.tag) ))
        elif nPars==6: model = ROOT.RooGenericPdf(pdfName, "(0.5*tanh((@0-@6)/@5) + .5)*pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1_%s"%self.tag), self.w.var("CMS_JJ_p2_%s"%self.tag), self.w.var("CMS_JJ_p3_%s"%self.tag), self.w.var("CMS_JJ_p4_%s"%self.tag), self.w.var("CMS_JJ_p5_%s"%self.tag), self.w.var("CMS_JJ_p6_%s"%self.tag)))

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))

    def addQCDShapeNoTag(self,name,variable,preconstrains,nPars=4):

        pdfName=name+"_"+self.tag
        
        MVV=variable
        if self.w.var(MVV) == None: self.w.factory(MVV+"[0,10000]")
        
        errs = [100,200,100,100,100]
        values = [9.28433e+00,1.03641e+01,2.35256e+00,4.17695e-01,1.00000e+01]
        f = ROOT.TFile.Open(preconstrains,'READ')
        parsG = [f.Get('p%i'%i) for i in range(1,nPars+1)]
        pars_val = [ROOT.Double(0.) for i in range(0,nPars)]       
        for i in range(1,nPars+1):
         x = ROOT.Double(0.)         
         pName="CMS_JJ_p%i"%i
         if nPars==6 and i==5:
            errUp=5000
            errDown=1
            #pars_val[i-1] = 1000
         elif nPars==6 and i==6:
            errUp=2500
            errDown=0
            #pars_val[i-1] = 1600
         else:
            errUp=errs[i-1]
            errDown=-errs[i-1]
            #pars_val[i-1] = values[i-1]
         parsG[i-1].GetPoint(0,x,pars_val[i-1])
         #errUp = pars_val[i-1]+parsG[i-1].GetErrorYhigh(0)*1000.
         #errDown = -pars_val[i-1]-parsG[i-1].GetErrorYlow(0)*1000.
         print i,pName,pars_val[i-1],parsG[i-1].GetErrorYhigh(0),parsG[i-1].GetErrorYlow(0),errUp,errDown
         self.w.factory("{name}[{val},{errDown},{errUp}]".format(name=pName,val=pars_val[i-1],errUp=errUp,errDown=errDown))
        
        if nPars==2: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2)", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2")))  
        elif nPars==3: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/pow(@0/13000., @2+@3*log(@0/13000.))", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3")))
        elif nPars==4: model = ROOT.RooGenericPdf(pdfName, "pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3"), self.w.var("CMS_JJ_p4")))
        elif nPars==5: model = ROOT.RooGenericPdf(pdfName, " pow(exp(-@0/13000.), @4) * pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3"), self.w.var("CMS_JJ_p4"), self.w.var("CMS_JJ_p5")))
        elif nPars==6: model = ROOT.RooGenericPdf(pdfName, "(0.5*tanh((@0-@6)/@5) + .5)*pow(1-@0/13000., @1)/ ( pow(@0/13000., @2+@3*log(@0/13000.)+@4*pow(log(@0/13000.),2)) )", ROOT.RooArgList(self.w.var(MVV), self.w.var("CMS_JJ_p1"), self.w.var("CMS_JJ_p2"), self.w.var("CMS_JJ_p3"), self.w.var("CMS_JJ_p4"), self.w.var("CMS_JJ_p5"), self.w.var("CMS_JJ_p6")))

        getattr(self.w,'import')(model,ROOT.RooFit.Rename(pdfName))
