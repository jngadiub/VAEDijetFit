import os,sys,time
from array import array

import ROOT
import CMS_lumi, tdrstyle

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
   
def plotPValue(xsec_scan):

	xmin = xsec_scan[0]
	xmax = xsec_scan[-1]+xsec_scan[-1]*0.1
    
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
	hrl_SM.GetXaxis().SetTitle("cross-section [pb]")

	x = array('d', xsec_scan)
	ys = array('d', [])
	yp = array('d',[])

	fin = open('results_accepted.txt')
	for l in fin.readlines():
 		l = l.split('\t')
 		yp.append(float(l[1]))
 		ys.append(float(l[2]))
	fin.close()
    
	nPoints=len(x)
	gp = ROOT.TGraph(nPoints,x,yp)
	gp.SetName("PValue")
	gp.SetLineColor(1)
	gp.SetMarkerColor(1)
	gp.SetMarkerStyle(20)
	gp.SetLineWidth(2)
	gp.SetMarkerSize(1.)

	ys2 = array('d', [])
	yp2 = array('d',[])

	fin = open('results_total.txt')
	for l in fin.readlines():
		l = l.split('\t')
 		yp2.append(float(l[1]))
 		ys2.append(float(l[2]))
	fin.close()
    
	gp2 = ROOT.TGraph(nPoints,x,yp2)
	gp2.SetName("PValue")
	gp2.SetLineColor(210)
	gp2.SetMarkerColor(210)
	gp2.SetMarkerStyle(20)
	gp2.SetLineWidth(2)
	gp2.SetMarkerSize(1.)

	ys3 = array('d', [])
	yp3 = array('d',[])

	fin = open('results_twoCatFit.txt')
	for l in fin.readlines():
		l = l.split('\t')
 		yp3.append(float(l[1]))
 		ys3.append(float(l[2]))
	fin.close()
    
	gp3 = ROOT.TGraph(nPoints,x,yp3)
	gp3.SetName("PValue")
	gp3.SetLineColor(ROOT.kPink)
	gp3.SetMarkerColor(ROOT.kPink)
	gp3.SetMarkerStyle(20)
	gp3.SetLineWidth(2)
	gp3.SetMarkerSize(1.)
	
	pvalues = [ ROOT.RooStats.SignificanceToPValue(i) for i in range(1,7) ]
	lines = [ ROOT.TF1("SLine_%i"%i,"%f"%pvalues[i-1],xmin,xmax) for i in range(1,7) ]
	for l in lines:
	 l.SetLineColor(ROOT.kRed)
	 l.SetLineWidth(2)
	 l.SetLineStyle(3)
	 
	bans = [ ROOT.TLatex(xmax*0.95,pvalues[i-1],("%i #sigma"%(i))) for i in range(1,7) ]
	for b in bans:
	 b.SetTextSize(0.028)
	 b.SetTextColor(2)

	legend = ROOT.TLegend(0.7510112,0.7183362,0.8502143,0.919833)
	legend.SetTextSize(0.032)
	legend.SetLineColor(0)
	legend.SetShadowColor(0)
	legend.SetLineStyle(1)
	legend.SetLineWidth(1)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetMargin(0.35)
	legend.AddEntry(gp,'Accepted (1-cat fit)','LP') 
	legend.AddEntry(gp3,'Accepted (2-cat fit)','LP') 
	legend.AddEntry(gp2,'Inclusive','LP') 

	gp.Draw('LP')
	gp2.Draw('LPsame')
	gp3.Draw('LPsame')
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
 
	canv.SaveAs("pvalue.png")
	time.sleep(1000)


if __name__ == "__main__":

 run = int(sys.argv[1])
 
 xsec = [i*0.0001 for i in range(0,15)]
 print xsec
 labels = ['total','accepted']
 mass = 3500.
 
 #if you have already run the scan, results are saved in txt files 
 if run == 0:
  plotPValue(xsec)
  sys.exit()

 x = array('d', xsec)
 ysig = array('d', [])
 ypvalue = array('d',[])
 ysig2 = array('d', [])
 ypvalue2 = array('d',[])

 #else run the scan for 1-category bump-hunt 
 for index in range(0,2):

	fname = 'results_{label}.txt'.format(label=labels[index])
	if os.path.exists(fname): os.system('rm %s'%fname)
	fout = open(fname,'w')
     
	for x in xsec:
 
 		cmd = 'python dijetfit.py --index {index} --xsec {xsec}'.format(index=index,xsec=x)
 		print "Executing:",cmd
 		os.system(cmd)
 		
 		tf = ROOT.TFile.Open('higgsCombinesignificance_{xsec}_{label}.Significance.mH{mass}.root'.format(xsec=x,label=labels[index],mass=int(mass)),'READ')
 		tree = tf.limit
 		tree.GetEntry(0) 		
 		ysig.append(tree.limit) 		
 		tf.Close()

 		tf = ROOT.TFile.Open('higgsCombinepvalue_{xsec}_{label}.Significance.mH{mass}.root'.format(xsec=x,label=labels[index],mass=int(mass)),'READ')
 		tree = tf.limit
 		tree.GetEntry(0) 		
 		ypvalue.append(tree.limit) 		
 		tf.Close()
 		
 		fout.write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x,pvalue=ypvalue[-1],sig=ysig[-1]))
 		 		 	
	fout.close()	
 
 #run scan for 2-categories bump hunt
 fname = 'results_twoCatFit.txt'
 if os.path.exists(fname): os.system('rm %s'%fname)
 fout = open(fname,'w')
 
 for x in xsec:
 
 	cmd = 'python dijetfit.py --index 2 --xsec {xsec} --twoCatFit'.format(xsec=x)
 	print "Executing:",cmd
 	os.system(cmd)

 	cmd = 'python dijetfit.py --index 1 --xsec {xsec} --twoCatFit'.format(xsec=x)
 	print "Executing:",cmd
 	os.system(cmd)

	tf = ROOT.TFile.Open('higgsCombinesignificance_{xsec}.Significance.mH{mass}.root'.format(xsec=x,mass=int(mass)),'READ')
	tree = tf.limit
	tree.GetEntry(0) 		
	ysig2.append(tree.limit) 		
	tf.Close()

	tf = ROOT.TFile.Open('higgsCombinepvalue_{xsec}.Significance.mH{mass}.root'.format(xsec=x,mass=int(mass)),'READ')
	tree = tf.limit
	tree.GetEntry(0) 		
	ypvalue2.append(tree.limit) 		
	tf.Close()
 		
	fout.write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=x,pvalue=ypvalue2[-1],sig=ysig2[-1]))
 		 	 			  
 print "******************** Results:"
 print 
 print " - xsec:"
 print
 print xsec
 print
 print " - significance (1-cat):"
 print
 print ysig
 print
 print " - pvalue (1-cat):"
 print
 print ypvalue
 print 
 print " - significance (2-cat):"
 print
 print ysig2
 print
 print " - pvalue (2-cat):"
 print
 print ypvalue2
 print 
 
 plotPValue(xsec)   
  
  
