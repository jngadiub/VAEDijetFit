import os
import ROOT
from array import array
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

from Utils import *


sample_names = {
                'GtoWW15naReco': r'$G(1.5 TeV)\to WW$ narrow',
                'GtoWW15brReco': r'$G(1.5 TeV)\to WW$ broad',
                'GtoWW25naReco': r'$G(2.5 TeV)\to WW$ narrow',
                'GtoWW25brReco': r'$G(2.5 TeV)\to WW$ broad',
                'GtoWW30naReco': r'$G(3.0 TeV)\to WW$ narrow',
                'GtoWW30brReco': r'$G(3.0 TeV)\to WW$ broad',
                'GtoWW35naReco': r'$G(3.5 TeV)\to WW$ narrow',
                'GtoWW35brReco': r'$G(3.5 TeV)\to WW$ broad',
                'GtoWW45naReco': r'$G(4.5 TeV)\to WW$ narrow',
                'GtoWW45brReco': r'$G(4.5 TeV)\to WW$ broad',
                }

sample_dirs = {
                'GtoWW15naReco': 'RSGraviton_WW_NARROW_13TeV_PU40_1.5TeV_reco',
                'GtoWW25naReco': 'RSGraviton_WW_NARROW_13TeV_PU40_2.5TeV_reco',
                'GtoWW30naReco': 'RSGraviton_WW_NARROW_13TeV_PU40_3.0TeV_reco',
                'GtoWW35naReco': 'RSGraviton_WW_NARROW_13TeV_PU40_3.5TeV_reco',
                'GtoWW45naReco': 'RSGraviton_WW_NARROW_13TeV_PU40_4.5TeV_reco',
                'GtoWW15brReco': 'RSGraviton_WW_BROAD_13TeV_PU40_1.5TeV_reco',
                'GtoWW25brReco': 'RSGraviton_WW_BROAD_13TeV_PU40_2.5TeV_reco',
                'GtoWW30brReco': 'RSGraviton_WW_BROAD_13TeV_PU40_3.0TeV_reco',
                'GtoWW35brReco': 'RSGraviton_WW_BROAD_13TeV_PU40_3.5TeV_reco',
                'GtoWW45brReco': 'RSGraviton_WW_BROAD_13TeV_PU40_4.5TeV_reco',
                }


def get_results_file(run_n, sig_id, sig_xsec, q, loss_id):
    sig_dir = make_dir_str(sample_dirs[sig_id], sig_xsec, run_n, loss_id)
    return os.path.join(sig_dir,'results_%s.txt'%q)


def make_significance_lines(xmin, xmax):
    # draw pvalue and significance horizontal lines
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

    return lines, bans


def plot_pvalue(run_n, sig_ids, sig_xsec, xsec_scan, quantiles, loss_id='', plot_name_suffix='', out_dir=''):

    xmin = xsec_scan[0]*1000.
    xmax = (xsec_scan[-1]+xsec_scan[-1]*0.1)*1000.
    
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
    palette = get_palette('mass')
    col = ROOT.TColor()

    # for each signal mass
    for (sig_n, sig_id) in enumerate(sig_ids):
        color_sig = col.GetColor(palette[sig_n])
        x = array('d', xsec_scan*1000.)
 
        # for each quantile in (total (=classic bump hunt), final=(AD bump hunt))
        for iq,q in enumerate(quantiles):

            line_style = ROOT.kSolid if q == 'final' else ROOT.kDashed
        
            ys = array('d', [])
            yp = array('d',[])

            fin = open(get_results_file(run_n, sig_id, sig_xsec, q, loss_id), 'r')
            for l in fin.readlines():
                l = l.split('\t')
                yp.append(float(l[1]))
                ys.append(float(l[2]))
            fin.close()
        
            nPoints=len(x)
            gp = ROOT.TGraph(nPoints,x,yp)
            gp.SetName("PValue_%s"%q)
            gp.SetLineColor(color_sig)
            gp.SetLineStyle(line_style)
            gp.SetMarkerColor(color_sig)
            gp.SetMarkerStyle(20)
            gp.SetLineWidth(2)
            gp.SetMarkerSize(1.)
            graphs.append(gp)

    lines, bans = make_significance_lines(xmin, xmax)

    legend = ROOT.TLegend(0.18,0.2183362,0.28,0.419833)
    legend.SetTextSize(0.032)
    legend.SetLineColor(0)
    legend.SetShadowColor(0)
    legend.SetLineStyle(1)
    legend.SetLineWidth(1)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetMargin(0.35)

    for g, sig_id in zip(range(0,len(graphs),2), sig_ids):
        legend.AddEntry(graphs[g], sample_names[sig_id], 'LP')

    graphs[0].Draw('LP')
    for g in range(1,len(graphs)): graphs[g].Draw("LPsame")
    for l in lines:
        print(l)
        l.Draw("same")

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
 
    canv.SaveAs(os.path.join(out_dir, "pvalue"+plot_name_suffix+'_loss_'+loss_id+".png"))


if __name__ == "__main__":

    run_n = 113
    resonance = 'br'
    sig_ids = ['GtoWW15'+resonance+'Reco', 'GtoWW25'+resonance+'Reco', 'GtoWW35'+resonance+'Reco', 'GtoWW45'+resonance+'Reco',]
    quantiles = ['total', 'final']
    sig_xsec = 100.
    xsec_scan = np.array([0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1])

    plot_pvalue_matplotlib(run_n, sig_ids, sig_xsec, xsec_scan, quantiles, plot_name_suffix='_all_resonances_'+resonance, out_dir='')
