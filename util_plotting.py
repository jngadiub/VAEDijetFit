from array import array
import random
import numpy as np
import os

import ROOT
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
import CMS_lumi, tdrstyle
tdrstyle.setTDRStyle()
ROOT.gROOT.SetBatch(True)
ROOT.RooRandom.randomGenerator().SetSeed(random.randint(0, 1e+6))


def get_palette(mode):
    
    palette = {}
    palette['gv'] = []

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']

    for c in colors:
        palette['gv'].append(c)
 
    return palette[mode]


def get_canvas(cname):

    tdrstyle.setTDRStyle()

    H_ref = 630
    W_ref = 600
    W = W_ref
    H  = H_ref

    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.04*W_ref

    canvas = ROOT.TCanvas(cname,cname,50,50,W,H)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillStyle(0)
    canvas.SetFrameBorderMode(0)
    canvas.SetLeftMargin( L/W+0.01 )
    canvas.SetRightMargin( R/W+0.03 )
    canvas.SetTopMargin( 0.07 ) #/T/H
    canvas.SetBottomMargin( B/H )
    #canvas.SetGrid()
    canvas.SetLogy()

    return canvas
    

def plotPValue(xsec_scan, quantiles, labels, plot_name_suffix='', out_dir=''):

    xsec_scan = np.asarray(xsec_scan, dtype=np.float32)

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
    # time.sleep(10)
