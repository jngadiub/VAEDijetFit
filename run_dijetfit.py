import os,sys,time
import numpy as np
from array import array
import optparse, json
import re

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

        print(out_dir)
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
 
    canv.SaveAs(os.path.join(out_dir, "pvalue"+plot_name_suffix+".png"))
    time.sleep(10)


if __name__ == "__main__":

    #python run_dijetfit.py --run -i /eos/user/k/kiwoznia/data/QR_results/events/qr_run_6160/env_run_0/poly_run_0/ --sig GtoWW35naReco.h5 --qcd qcdSigAll.h5 --xsec 100 -M 3500 --res na -C

    parser = optparse.OptionParser()
    parser.add_option("--run","--run", dest="run", default=False, action="store_true", help="Run scan")
    parser.add_option("-n","-n",dest="run_n", type=int, default=0, help="Experiment number")
    parser.add_option("-M","-M", dest="mass", type=float, default=3500., help="Injected signal mass")
    parser.add_option("-i","--inputDir",dest="inputDir", default='./', help="directory with all quantiles h5 files")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcd.h5', help="QCD h5 file")
    parser.add_option("--sig","--sig", dest="sigFile", default='signal.h5', help="Signal h5 file")
    parser.add_option("-xsec", "--sigxsec", dest="sigXsec", default=100, help="true signal cross-section")
    parser.add_option("--res", "--res", dest="sigRes", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
    parser.add_option('-C', dest="correlateB",action="store_true",help="Coorelate background shape among quantiles")
    (options,args) = parser.parse_args()

    run = options.run
    mass = options.mass
    sigFile = options.sigFile
    qcdFile = options.qcdFile
    inputDir = options.inputDir
    sigRes = options.sigRes
    qr_run = re.findall(r'qr_run_(\d+)',options.inputDir)[0]
    xsec = np.array(get_xsec_scan_from_injection(options.sigXsec)) # pb

    # distinctive run string
    run_str = make_run_str(sig_name=options.sigFile, sig_xsec=options.sigXsec, run_n=qr_run)
    out_dir = run_str[1:]
    os.system('mkdir %s'%out_dir)
    print(run_str[1:])

    if len(xsec) == 0:
        print "ERROR: set the cross sections to scan for signal",sigFile,"in the files_count.json file!"
        sys.exit()

    quantiles = ['q05', 'q10', 'q30', 'q50', 'q70', 'q100', 'total']
    labels = ['q:0-5%','q:5-10%','q:10-30%','q:30-50%','q:50-70%','q:70-100%','bump hunt']

    #if you have already run the scan, results are saved in txt files 
    if run == 0:
        plotPValue(xsec, quantiles + ['final'], labels + ['AD bump hunt'], run_str, out_dir=out_dir)
        print("NOT CHECK OUTPUT FOLDER",out_dir)
        sys.exit()

    #now run the scan
    x = array('d', xsec)
    ysig = {}
    ypvalue = {}
    outfiles = []
    for q in quantiles:
        ysig[q] = []
        ypvalue[q] = []
        outfiles.append(open(out_dir+'/results_%s.txt'%q,'w'))

    for x in xsec:

        cmd = "python dijetfit.py -i {inputdir} --sig {sigfile} --qcd {qcdfile} --xsec {xsec} -M {mass} --res {res} --out {out_dir}".format(inputdir=inputDir, xsec=x, sigfile=sigFile, qcdfile=qcdFile, mass=mass, res=sigRes, out_dir=out_dir)
        if x!=0: cmd+=' -l' #assuming first xsec is zero so after that you can load input data with option -l
        if options.correlateB == True: cmd += ' -C'
        print cmd
        os.system(cmd)

        for iq,q in enumerate(quantiles):

            cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n significance_{xsec}_{label}'.format(out_dir=out_dir, xsec=x, label=q, mass=int(mass))
            print cmd
            os.system(cmd)

            cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{label}.root -m {mass} -n pvalue_{xsec}_{label} --pvalue'.format(out_dir=out_dir, xsec=x, label=q, mass=int(mass))
            print cmd
            os.system(cmd)
                
            tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}_{label}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=x,label=q,mass=int(mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ysig[q].append(tree.limit)
            print "Xsec",x,"quantile",q,"significance",ysig[q][-1]       
            tf.Close()

            tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}_{label}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=x,label=q,mass=int(mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ypvalue[q].append(tree.limit)        
            tf.Close()

            outfiles[iq].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x,pvalue=ypvalue[q][-1],sig=ysig[q][-1]))
                 
    for iq,q in enumerate(quantiles): outfiles[iq].close() 
 
    ysig['combo'] = []
    ypvalue['combo'] = []
    outfiles.append(open(out_dir+'/results_final.txt','w'))

    for x in xsec:

        cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n significance_{xsec}'.format(out_dir=out_dir, xsec=x,label='final',mass=int(mass))
        print cmd
        os.system(cmd)

        cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n pvalue_{xsec} --pvalue'.format(out_dir=out_dir, xsec=x,label='final',mass=int(mass))
        print cmd
        os.system(cmd)
            
        tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=x,mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ysig['combo'].append(tree.limit)             
        print "Xsec",x,"COMBO significance",ysig['combo'][-1]        
        tf.Close()

        tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=x,mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ypvalue['combo'].append(tree.limit)          
        tf.Close()

        outfiles[-1].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x, pvalue=ypvalue['combo'][-1], sig=ysig['combo'][-1]))  
 
    outfiles[-1].close()
    
    print ysig
    print ypvalue
   
    plotPValue(xsec, quantiles + ['final'], labels + ['AD bump hunt'], run_str, out_dir=out_dir)
    print("NOT CHECK OUTPUT FOLDER",out_dir)
  
