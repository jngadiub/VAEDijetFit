
import ROOT
import tdrstyle
import json


# define color palett for pvalue plotting
palette = {}
palette['gv'] = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
                '#bcbd22', '#17becf']
palette['mass'] = ['#3E96A1', '#EC4E20', '#FF9505', '#713E5A']


def get_palette(mode):
 return palette[mode]

 
def getBinning(binsMVV,minx,maxx,bins):
    l=[]
    if binsMVV=="":
        print(" do this")
        print(binsMVV)
        for i in range(0,bins+1):
            l.append(minx + i* (maxx - minx)/bins)
    else:
        print("dot that")
        print(binsMVV)
        s = binsMVV.split(",")
        for w in s:
            l.append(int(w))
    return l

def truncate(binning,mmin,mmax):
    res=[]
    for b in binning:
        if b >= mmin and b <= mmax:
            res.append(b)
    return res

def make_run_str(sig_name, sig_xsec=10, run_n=0, loss_id='rk5'):
    return '_' + sig_name[:-3] + '_xsec' + str(sig_xsec) + '_run' + str(run_n) + '_loss_' + loss_id 


def make_dir_str(sig_name, sig_xsec=10, run_n=0, loss_id='rk5'):
    return sig_name + '_xsec' + str(sig_xsec) + '_run' + str(run_n) + '_loss_' + loss_id 


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

def get_xsec_scan(filename):

    with open('files_count.json') as f:
        data = json.load(f)

    for k in data.keys():
        if k in filename or k.replace('_EXT','') in filename: return data[k][2]
   
