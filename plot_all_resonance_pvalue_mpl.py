# plots pvalue plot with matplotlib,
# need separate script because ROOT can only be used with sourced environment and python 2
# and mpl hep style has a bug in python 2

import os
from array import array
import numpy as np
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep

# define color palett for pvalue plotting
palette = {}
palette['gv'] = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                '#bcbd22', '#17becf']
palette['mass'] = ['#3E96A1', '#EC4E20', '#FF9505', '#713E5A']


def get_palette(mode):
 return palette[mode]

resonance_names = {
    'na' : r'Narrow $G\to WW$',
    'br' : r'Broad $G\to WW$',
}


sample_mass_names = {
                'GtoWW15naReco': r'$1.5 \, TeV$',
                'GtoWW15brReco': r'$1.5 \, TeV$',
                'GtoWW25naReco': r'$2.5 \, TeV$',
                'GtoWW25brReco': r'$2.5 \, TeV$',
                'GtoWW30naReco': r'$3.0 \, TeV$',
                'GtoWW30brReco': r'$3.0 \, TeV$',
                'GtoWW35naReco': r'$3.5 \, TeV$',
                'GtoWW35brReco': r'$3.5 \, TeV$',
                'GtoWW45naReco': r'$4.5 \, TeV$',
                'GtoWW45brReco': r'$4.5 \, TeV$',
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

def make_dir_str(sig_name, sig_xsec=10, run_n=0):
    return sig_name + '_xsec' + str(sig_xsec) + '_run' + str(run_n) 


def get_results_file(run_n, sig_id, sig_xsec, q):
    sig_dir = make_dir_str(sample_dirs[sig_id], sig_xsec, run_n)
    return os.path.join(sig_dir,'results_%s.txt'%q)


def compute_pval_from_zscore(zscores):
    return 1 - stats.norm.cdf(zscores)


def plot_pvalue_matplotlib(run_n, resonance, sig_ids, sig_xsec, xsec_scan, quantiles, loss_id, plot_name_suffix, fig_dir):
    
    # Load CMS style sheet
    plt.style.use(hep.style.CMS)

    xmin = xsec_scan[0]*1000.
    xmax = (xsec_scan[-1]+xsec_scan[-1]*0.1)*1000

    palette = get_palette('mass')
    fig = plt.figure() # figsize=(5, 5)

    # for each signal mass
    for (sig_n, sig_id) in enumerate(sig_ids):

        color_sig = palette[sig_n]
        x = array('d', xsec_scan*1000.)

        # for each quantile in (total (=classic bump hunt), final=(AD bump hunt))
        for iq,q in enumerate(quantiles):

            line_style = 'solid' if q == 'final' else 'dashed'
        
            ys = array('d', [])
            yp = array('d', [])

            fin = open(get_results_file(run_n, sig_id, sig_xsec, q), 'r')
            for l in fin.readlines():
                l = l.split('\t')
                yp.append(float(l[1]))
                ys.append(float(l[2]))
            fin.close()

            plt.semilogy(x, yp, color=color_sig, linestyle=line_style, marker='.')

    # add pvalue & significance lines
    pvals = compute_pval_from_zscore(np.arange(1,7))
    for s, p in enumerate(pvals):
        plt.axhline(y=p, color='gray', linestyle='dotted', linewidth=0.8)
        plt.text(xmax*0.9, p, r'$'+str(s+1)+r' \sigma$', color='gray')

    plt.xlabel('cross-section [fb]', loc='right')
    plt.ylabel('p-value', loc='top')
    plt.ylim(bottom=1e-12, top=1)
    plt.xlim(left=0)

    # set yticks manually
    plt.gca().tick_params(direction='in', which='both')
    plt.gca().minorticks_on()
    locmaj = matplotlib.ticker.LogLocator(base=10, numticks=15) 
    plt.gca().yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=10.0, subs=(0.2,0.4,0.6,0.8), numticks=15)
    plt.gca().yaxis.set_minor_locator(locmin)
    plt.gca().yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

    # plt.legend(loc='best', frameon=False)
    lines = plt.gca().get_lines()
    legend1 = plt.legend([lines[i] for i in [0,2,4,6]], [sample_mass_names[sig_id] for sig_id in sig_ids], loc='best', frameon=False, markerscale=3, handlelength=1, title=resonance_names[resonance])
    legend2 = plt.legend([lines[i] for i in [0,1]], ['bump hunt', 'AD bump hunt'], loc='center left', frameon=False, markerscale=0)
    for leg in legend2.legendHandles: 
        leg.set_color('black')
    plt.gca().add_artist(legend1)
    plt.gca().add_artist(legend2)

    print('writing ROC plot to {}'.format(fig_dir))
    fig.savefig(os.path.join(fig_dir, "pvalue"+ plot_name_suffix + '_loss_'+loss_id + '_mpl.png'), bbox_inches='tight')
    plt.close(fig)



if __name__ == "__main__":

    run_n = 113
    resonance = 'na'
    sig_ids = ['GtoWW15'+resonance+'Reco', 'GtoWW25'+resonance+'Reco', 'GtoWW35'+resonance+'Reco', 'GtoWW45'+resonance+'Reco',]
    quantiles = ['total', 'final']
    sig_xsec = 100.
    xsec_scan = np.array([0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1])

    plot_pvalue_matplotlib(run_n, resonance, sig_ids, sig_xsec, xsec_scan, quantiles, loss_id='rk5_05', plot_name_suffix='_all_resonances_'+resonance, fig_dir='')
