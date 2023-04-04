import os,sys,time
import numpy as np
from array import array
import optparse, json

import ROOT
import CMS_lumi, tdrstyle

from Utils import *


def plotPValue(xsec_scan, quantiles, labels, plot_name_suffix='', out_dir=''):

    scale = 1000. #injected xsec is pb but we want to plot in fb
    xmin = xsec_scan[0]*scale
    xmax = (xsec_scan[-1]+xsec_scan[-1]*0.1)*scale
    
    canv = get_canvas("c_Significance")
    canv.cd()

    hrl_SM = canv.DrawFrame(xmin,1e-12,xmax, 1) 
    ytitle = "p-value"
    hrl_SM.GetYaxis().SetTitle(ytitle)
    hrl_SM.GetYaxis().SetTitleSize(0.045)
    hrl_SM.GetXaxis().SetTitleSize(0.045)
    hrl_SM.GetXaxis().SetLabelSize(0.03)
    hrl_SM.GetYaxis().SetLabelSize(0.03)
    hrl_SM.GetYaxis().SetTitleOffset(1.2)
    hrl_SM.GetXaxis().SetTitleOffset(1.1)
    hrl_SM.GetXaxis().SetTitle("cross-section [fb]")

    graphs = []
    palette = get_palette('gv')
    col = ROOT.TColor()
 
    for iq,q in enumerate(quantiles):
    
        x = array('d', xsec_scan*scale)
        ys = array('d', [])
        yp = array('d',[])

        fin = open(os.path.join(out_dir,'results_%s.txt'%q), 'r')
        for l in fin.readlines():
            l = l.split('\t')
            yp.append(float(l[1]))
            ys.append(float(l[2]))
        fin.close()
    
        print(iq,q,ys)
        print(iq,q,yp)
        nPoints=len(x)
        gp = ROOT.TGraph(nPoints,x,yp)
        gp.SetName("PValue_%s"%q)
        if q!='total':
            gp.SetLineColor(col.GetColor(palette[iq]))
            gp.SetMarkerColor(col.GetColor(palette[iq]))
        else: 
            gp.SetLineColor(1)
            gp.SetMarkerColor(1)
        gp.SetMarkerStyle(20)
        gp.SetLineWidth(2)
        gp.SetMarkerSize(1.)
        graphs.append(gp)

    #JEN -- to add the q90 only category on top
    #fin = open(os.path.join("eos_vae/scans_XToYYprimeTo4Q_MX2000_MY400_MYprime170/scan_run1_signal_XToYYprimeTo4Q_MX2000_MY400_MYprime170_narrowReco",'results_q90.txt'), 'r')
    #x = array('d', xsec_scan*scale)
    #ys = array('d', [])
    #yp = array('d',[])
    #for l in fin.readlines():
    #        l = l.split('\t')
    #        yp.append(float(l[1]))
    #        ys.append(float(l[2]))
    #fin.close()
    
    #print("q90-only",ys)
    #print("q90-only",yp)
    #nPoints=len(x)
    #gp = ROOT.TGraph(nPoints,x,yp)
    #gp.SetName("PValue_q90_only")
    #gp.SetLineColor(col.GetColor(palette[len(quantiles)+1]))
    #gp.SetMarkerColor(col.GetColor(palette[len(quantiles)+1]))
    #gp.SetMarkerStyle(20)
    #gp.SetLineWidth(2)
    #gp.SetMarkerSize(1.)
    #graphs.append(gp)
    #JEN

    pvalues = [ ROOT.RooStats.SignificanceToPValue(i) for i in range(1,7) ]
    lines = [ ROOT.TLine(xmin,pvalues[i-1],xmax,pvalues[i-1]) for i in range(1,7) ]
    for l in lines:
        l.SetLineColor(ROOT.kGray+1)
        l.SetLineWidth(2)
        l.SetLineStyle(ROOT.kDashed)
     
    bans = [ ROOT.TLatex(xmax*0.93,pvalues[i-1],("%i #sigma"%(i))) for i in range(1,7) ]
    for b in bans:
        b.SetTextSize(0.028)
        b.SetTextColor(ROOT.kGray+1)

    legend = ROOT.TLegend(0.18,0.2183362,0.28,0.419833)
    legend.SetTextSize(0.032)
    legend.SetLineColor(0)
    legend.SetShadowColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetMargin(0.35)

    for iq,q in enumerate(quantiles): legend.AddEntry(graphs[iq],labels[iq],'LP')
    #legend.AddEntry(graphs[iq+1],labels[iq+1],'LP') #JEN -- to add the q90 only category on top

    graphs[0].Draw('LP')
    for g in range(1,len(graphs)): graphs[g].Draw("LPsame")
    for l in lines: l.Draw("same")

    for b in bans: b.Draw()

    legend.Draw()
    canv.Update() 
    canv.cd()
    canv.Update()
    canv.RedrawAxis()
    canv.RedrawAxis("g")
    frame = canv.GetFrame()
    frame.Draw() 
    canv.cd()
    canv.Update()
 
    canv.SaveAs(os.path.join(out_dir, "pvalue_"+plot_name_suffix+".png"))
    time.sleep(10)


if __name__ == "__main__":

    #python run_dijetfit.py --i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 --qcd bkg.h5 -M 2000 -n 1 --run --config 4 -C

    parser = optparse.OptionParser()
    parser.add_option("--run","--run", dest="run", default=False, action="store_true", help="Run scan")
    parser.add_option("-n","-n",dest="run_n", type=int, default=0, help="Experiment number")
    parser.add_option("-M","-M", dest="mass", type=float, default=3500., help="Injected signal mass")
    parser.add_option("-i","--inputDir",dest="inputDir", default='./', help="directory with all quantiles h5 files")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcd.h5', help="QCD h5 file")
    parser.add_option("--sig","--sig",dest="signal",type=str,default='XToYYprimeTo4Q_MX2000_MY80_MYprime170',help="Signal name")
    parser.add_option("-C","--correlate",dest="correlateB",action="store_true",help="Coorelate background shape among quantiles")
    parser.add_option("--config","--config", dest="config", type=int, default=4, help="quantiles config")
    parser.add_option("--init","--init",dest="init",action="store_true",help="True to manually inject signal in spectrum and retrieve observed significance; False to let combine do the emulatio and retrieve expected significance")
    (options,args) = parser.parse_args()

    init = options.init
    run = options.run
    mass = options.mass
    signal = options.signal
    qcdFile = options.qcdFile
    inputDir = options.inputDir
    xsec = np.array(get_xsec_scan(options.signal)) # pb

    out_dir = "scan_run"+str(options.run_n)+"_"+signal
    os.system('mkdir %s'%out_dir)

    if len(xsec) == 0:
        print "ERROR: set the cross sections to scan for signal",sigFile,"in the files_count.json file!"
        sys.exit()

    if options.config==1:
        quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'total']
        labels = ['q:70-100%', 'q:50-70%', 'q:30-50%', 'q:10-30%', 'q:0-10%', 'bump hunt']
    elif options.config==2:    
        quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'q95', 'total']
        labels = ['q:70-100%', 'q:50-70%', 'q:30-50%', 'q:10-30%', 'q:5-10%', 'q:0-5%', 'bump hunt']
    elif options.config==3:    
        quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'q95', 'q99', 'total']
        labels = ['q:70-100%', 'q:50-70%', 'q:30-50%', 'q:10-30%', 'q:5-10%', 'q:1-5%', 'q:0-1%', 'bump hunt']
    elif options.config==4:    
        quantiles = ['q90', 'q95', 'q99', 'total']
        #labels = ['q: 90 - 95%', 'Q: 95 - 99%', 'Q: > 99%', 'Inclusive Fit', "Q: > 90%"]
        labels = ['q: 90 - 95%', 'Q: 95 - 99%', 'Q: > 99%', 'Inclusive Fit']
    elif options.config==5:    
        quantiles = ['q95', 'q99', 'total']
        labels = ['q:1-5%', 'q:0-1%', 'bump hunt']
    elif options.config==6:    
        quantiles = ['q70', 'q90', 'q95', 'q99', 'total']
        labels = ['q:10-30%', 'q:5-10%', 'q:1-5%', 'q:0-1%', 'bump hunt']
    elif options.config==7:    
        quantiles = ['q50', 'q70', 'q90', 'q95', 'q99', 'total']
        labels = ['q:30-50%', 'q:10-30%', 'q:5-10%', 'q:1-5%', 'q:0-1%', 'bump hunt']

    #if you have already run the scan, results are saved in txt files 
    if run == 0:
        if init == True: out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
        plotPValue(xsec, quantiles + ['final'], labels + ['AD Combined Fit'], signal, out_dir=out_dir)
        print("NOT CHECK OUTPUT FOLDER",out_dir)
        sys.exit()

    if init == True:
        #prepare xsec=0 workspace
        cmd = "python dijetfit.py --xsec 0.0 -M {mass} -i {inputdir} --sig signal_{signal}_narrowReco.h5 --qcd {qcdfile} --out scan_run{run}_{signal}_xsec0_init --config {config} --run_toys".format(run=options.run_n, inputdir=inputDir, signal=signal, qcdfile=qcdFile, mass=mass, config=options.config)
        if options.correlateB == True: cmd += ' -C'
        print(cmd)
        os.system(cmd)

        #prepare xsec>0 workspace [this is needed so that the --expectSignal in combine is equal to the cross section]
        cmd = "python dijetfit.py --xsec 1.0 -M {mass} -i {inputdir} --sig signal_{signal}_narrowReco.h5 --qcd {qcdfile} --out scan_run{run}_{signal}_xsec1_init --config {config} --run_toys".format(run=options.run_n, inputdir=inputDir, signal=signal, qcdfile=qcdFile, mass=mass, config=options.config)
        if options.correlateB == True: cmd += ' -C'
        print(cmd)
        os.system(cmd)

    #now run the scan
    x = array('d', xsec)
    ysig = {}
    ypvalue = {}
    outfiles = []
    for q in quantiles:
        ysig[q] = []
        ypvalue[q] = []
        if init: out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
        outfiles.append(open(out_dir+'/results_%s.txt'%q,'w'))

    for x in xsec:

        qcdfile = qcdFile
        exp_sig = 0
        xsec_init = 0.0
        if init == True:
            exp_sig = x
            out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
        if x>0: 
            exp_sig = 1
            xsec_init = 1.0
            qcdfile = qcdFile.replace('.h5',"_"+signal+"_narrowReco_"+str(x)+'.h5') #take background file with signal injected in the QR
            if init == True:
                exp_sig = x
                out_dir = "scan_run{run}_{signal}_xsec1_init".format(run=options.run_n, signal=signal)
                qcdfile = qcdFile
        
        if init == False: #manual signal injection in background spectrum
            cmd = "python dijetfit.py -i {inputdir} --sig signal_{signal}_narrowReco.h5 --qcd {qcdfile} --xsec {xsec} -M {mass} --out {out_dir}/xsec{xsec} --config {config}".format(inputdir=inputDir, xsec=x, signal=signal, qcdfile=qcdfile, mass=mass, out_dir=out_dir,config=options.config)
            if options.correlateB == True: cmd += ' -C' 
            print cmd
            os.system(cmd)
        
        for iq,q in enumerate(quantiles):

            cmd = 'cd {out_dir}/xsec{XSEC} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{xsec}_{label}'.format(out_dir=out_dir, XSEC=x, xsec=x*1000, label=q, mass=int(mass))
            if init == True: cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{expsig}_{label} --toysFreq -t -1 --expectSignal={expsig}'.format(expsig=exp_sig,out_dir=out_dir, xsec=xsec_init*1000, label=q, mass=int(mass))
            print cmd
            os.system(cmd)

            cmd = 'cd {out_dir}/xsec{XSEC} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n pvalue_{xsec}_{label} --pvalue'.format(out_dir=out_dir, XSEC=x, xsec=x*1000, label=q, mass=int(mass))
            if init == True: cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n pvalue_{expsig}_{label} --pvalue --toysFreq -t -1 --expectSignal={expsig}'.format(expsig=exp_sig, out_dir=out_dir, xsec=xsec_init*1000, label=q, mass=int(mass))
            print cmd
            os.system(cmd)
                
            if init == False: tf = ROOT.TFile.Open('{out_dir}/xsec{XSEC}/higgsCombinesignificance_{xsec}_{label}.Significance.mH{mass}.root'.format(out_dir=out_dir, XSEC=x, xsec=x*1000,label=q,mass=int(mass)),'READ')
            else: tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{expsig}_{label}.Significance.mH{mass}.root'.format(expsig=exp_sig, out_dir=out_dir,label=q,mass=int(mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ysig[q].append(tree.limit)
            print "Xsec",x,"quantile",q,"significance",ysig[q][-1]       
            tf.Close()

            if init == False: tf = ROOT.TFile.Open('{out_dir}/xsec{XSEC}/higgsCombinepvalue_{xsec}_{label}.Significance.mH{mass}.root'.format(out_dir=out_dir, XSEC=x, xsec=x*1000,label=q,mass=int(mass)),'READ')
            else: tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{expsig}_{label}.Significance.mH{mass}.root'.format(expsig=exp_sig, out_dir=out_dir, label=q,mass=int(mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ypvalue[q].append(tree.limit)        
            tf.Close()

            outfiles[iq].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x,pvalue=ypvalue[q][-1],sig=ysig[q][-1]))
                 
    for iq,q in enumerate(quantiles): outfiles[iq].close() 
 
    ysig['combo'] = []
    ypvalue['combo'] = []
    if init: out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
    outfiles.append(open(out_dir+'/results_final.txt','w'))

    for x in xsec:

        exp_sig = 0
        xsec_init = 0.0
        if init == True:
            exp_sig = x
            out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
        if x>0: 
            exp_sig = 1
            xsec_init = 1.0
            if init == True:
                exp_sig = x
                out_dir = "scan_run{run}_{signal}_xsec1_init".format(run=options.run_n, signal=signal)
            
        if init == False: cmd = 'cd {out_dir}/xsec{XSEC} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{xsec}'.format(out_dir=out_dir, XSEC=x, xsec=x*1000.,label='final',mass=int(mass))
        else: cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{expsig} --toysFreq -t -1 --expectSignal={expsig}'.format(expsig=exp_sig, xsec=xsec_init*1000, out_dir=out_dir, label='final',mass=int(mass))
        print cmd
        os.system(cmd)

        if init == False: cmd = 'cd {out_dir}/xsec{XSEC} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n pvalue_{xsec} --pvalue'.format(out_dir=out_dir, XSEC=x, xsec=x*1000.,label='final',mass=int(mass))
        else: cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n pvalue_{expsig} --pvalue  --toysFreq -t -1 --expectSignal={expsig}'.format(expsig=exp_sig, xsec=xsec_init*1000, out_dir=out_dir, label='final',mass=int(mass))
        print cmd
        os.system(cmd)
            
        if init == False: tf = ROOT.TFile.Open('{out_dir}/xsec{XSEC}/higgsCombinesignificance_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, XSEC=x, xsec=x*1000.,mass=int(mass)),'READ')
        else: tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{expsig}.Significance.mH{mass}.root'.format(expsig=exp_sig, out_dir=out_dir, mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ysig['combo'].append(tree.limit)             
        print "Xsec",x,"COMBO significance",ysig['combo'][-1]        
        tf.Close()

        if init == False: tf = ROOT.TFile.Open('{out_dir}/xsec{XSEC}/higgsCombinepvalue_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, XSEC=x, xsec=x*1000.,mass=int(mass)),'READ')
        else: tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{expsig}.Significance.mH{mass}.root'.format(expsig=exp_sig, out_dir=out_dir, mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ypvalue['combo'].append(tree.limit)          
        tf.Close()

        outfiles[-1].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x, pvalue=ypvalue['combo'][-1], sig=ysig['combo'][-1]))  
 
    outfiles[-1].close()
    
    print ysig
    print ypvalue
   
    if init: out_dir = "scan_run{run}_{signal}_xsec0_init".format(run=options.run_n, signal=signal)
    plotPValue(xsec, quantiles + ['final'], labels + ['AD bump hunt'], signal, out_dir=out_dir)
    print("NOT CHECK OUTPUT FOLDER",out_dir)
  
