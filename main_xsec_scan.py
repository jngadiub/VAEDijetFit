import numpy as np
import optparse
import json
import pathlib2 as pathli  # python 2 backport
import os
import ROOT
from collections import defaultdict
from util_plotting import *


if __name__ == "__main__":

    #python main_xsec_scan.py -M 3500 --sig GtoWW35naReco.h5 --qcd qcdSigAll.h5 --res na -C --plt

    parser = optparse.OptionParser()
    parser.add_option("-M","-M", dest="mass", type=float, default=3500., help="Injected signal mass")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcdSigAll.h5', help="QCD signalregion h5 file")
    parser.add_option("--sig","--sig", dest="sigFile", default='signal.h5', help="Signal h5 file")
    parser.add_option("--res", "--res", dest="sigRes", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
    parser.add_option('-C', dest="correlateB",action="store_true",help="Coorelate background shape among quantiles")
    parser.add_option('--signi', dest='signi',action='store_true', help='run only the significance tests')
    parser.add_option('--plt', dest='plt_only',action='store_true', help='plot only xsec scan results')
    parser.add_option('--gof', dest='do_gof',action='store_true', help='run the goodness of fit test')
    (opts,args) = parser.parse_args()


    dict_g35na = {
        0  : 416,
        50 : 6161,
        100: 6160
    }
    experiment_dict = { 'GtoWW35naReco.h5': dict_g35na }
    dd = experiment_dict[opts.sigFile]

    quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'total']
    labels = ['q:0-30%','q:30-50%','q:50-70%','q:70-90%','q:90-100%','bump hunt']

    # base output dir for signal
    out_dir_base = opts.sigFile[:-len('Reco.h5')]
    pathli.Path(out_dir_base).mkdir(parents=True, exist_ok=True)

    if not opts.plt_only:

        # open the outfiles, ysig and ypvalue datastructs to collect results across xsecs for each quantile
        outfiles = {q: open('{}/results_{}.txt'.format(out_dir_base,q),'w') for q in quantiles}
        ysig = defaultdict(list)
        ypvalue = defaultdict(list)

        for xsec, qr_run in dd.items():

            # ****************************************************
            #                       BUMP HUNT
            # ****************************************************

            ## paths

            # distinct input dir for each cross section
            in_dir = '/eos/user/k/kiwoznia/data/QR_results/events/qr_run_'+str(dd[xsec])+'/env_run_0/poly_run_0/'
            # distinct output subdir for each cross section
            out_dir = os.path.join(out_dir_base,'xsec'+str(int(xsec))+'_qr'+str(int(qr_run)))
            if opts.correlateB is False:
                out_dir += '_noco'
            xsec = float(xsec)

            ## run dijetfit
            if not opts.signi:
                cmd = "python dijetfit.py -i {in_dir} --sig {sigfile} --qcd {qcdfile} --xsec {xsec} -M {mass} --res {res} --out {out_dir}".format(in_dir=in_dir, xsec=xsec, sigfile=opts.sigFile, qcdfile=opts.qcdFile, mass=opts.mass, res=opts.sigRes, out_dir=out_dir)
                if options.correlateB is True: 
                    cmd += ' -C'
                print(cmd)
                os.system(cmd)

            ## calculate significance
            for q in quantiles:

                cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{quantile}.root -m {mass} -n significance_{xsec}_{quantile}'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(opts.mass))
                print(cmd)
                os.system(cmd)

                cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{quantile}.root -m {mass} -n pvalue_{xsec}_{quantile} --pvalue'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(opts.mass))
                print(cmd)
                os.system(cmd)
                    
                tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}_{quantile}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(opts.mass)), 'READ')
                tree = tf.limit
                tree.GetEntry(0)         
                ysig[q].append(tree.limit)
                print("Xsec {}, quantile {}, significance {}".format(xsec,q,ysig[q][-1]))       
                tf.Close()

                tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}_{quantile}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(opts.mass)), 'READ')
                tree = tf.limit
                tree.GetEntry(0)         
                ypvalue[q].append(tree.limit)        
                tf.Close()

                outfiles[q].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=xsec,pvalue=ypvalue[q][-1],sig=ysig[q][-1]))

            # end for each quantile

        # end for each xsec
        for q in quantiles: outfiles[q].close()

        ### combination

        outfiles['combo'] = open(out_dir_base+'/results_final.txt','w') # save final quantile-combined xsec scan results to signal root directory

        for xsec, qr_run in dd.items():

            out_dir = os.path.join(out_dir_base,'xsec'+str(int(xsec))+'_qr'+str(int(qr_run)))
            xsec = float(xsec)

            cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n significance_{xsec}'.format(out_dir=out_dir, xsec=xsec,label='final',mass=int(opts.mass))
            print cmd
            os.system(cmd)

            cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n pvalue_{xsec} --pvalue'.format(out_dir=out_dir, xsec=xsec,label='final',mass=int(opts.mass))
            print cmd
            os.system(cmd)
                
            tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec,mass=int(opts.mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)             
            ysig['combo'].append(tree.limit)             
            print "Xsec",xsec,"COMBO significance",ysig['combo'][-1]        
            tf.Close()

            tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec,mass=int(opts.mass)),'READ')
            tree = tf.limit
            tree.GetEntry(0)             
            ypvalue['combo'].append(tree.limit)          
            tf.Close()

            outfiles['combo'].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=xsec, pvalue=ypvalue['combo'][-1], sig=ysig['combo'][-1]))  
     
        outfiles['combo'].close()
        
        print ysig
        print ypvalue

    # end if plot only

    plotPValue(list(dd.keys()), quantiles + ['final'], labels + ['AD bump hunt'], '_'.join([str(qr) for qr in dd.values()]), out_dir=out_dir_base)
    print("CHECK OUTPUT FOLDER",out_dir_base)

    # ****************************************************
    #                       GOF
    # ****************************************************

    if opts.do_gof is True:

        # e.g. python gof.py -i GtoWW35na/xsec100_qr6160 -N 1000 -T -C
        
        toys_n = 1000

        for xsec, qr_run in dd.items():

            dijet_out_dir = os.path.join(out_dir_base,'xsec'+str(int(xsec))+'_qr'+str(int(qr_run)))
            if opts.correlateB is False:
                    dijet_out_dir += '_noco'

            cmd = 'python gof.py -i {} --xsec {} -N {} -T'.format(dijet_out_dir, xsec, toys_n)
            if options.correlateB == True: 
                cmd += ' -C'

            print(cmd)
            os.system(cmd)

            # open the file produced by gof and append to xsec results collection
            

