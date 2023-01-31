import numpy as np
import optparse
import json
import pathlib
import os


if __name__ == "__main__":

    #python main_xsec_scan.py -M 3500 --sig GtoWW35naReco.h5 --qcd qcdSigAll.h5 --res na -C --signi

    parser = optparse.OptionParser()
    parser.add_option("-M","-M", dest="mass", type=float, default=3500., help="Injected signal mass")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcdSigAll.h5', help="QCD signalregion h5 file")
    parser.add_option("--sig","--sig", dest="sigFile", default='signal.h5', help="Signal h5 file")
    parser.add_option("--res", "--res", dest="sigRes", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
    parser.add_option('-C', dest="correlateB",action="store_true",help="Coorelate background shape among quantiles")
    parser.add_option('--signi', dest='signi',action='store_true', help='run only the significance tests')
    (opts,args) = parser.parse_args()


    dict_g35na = {
        0  : 416,
        50 : 6161,
        100: 6160
    }
    experiment_dict = { 'GtoWW35naReco.h5': dict_g35na }
    dd = experiment_dict[opt.sigFile]

    quantiles = ['q0', 'q30', 'q50', 'q70', 'q90', 'total']

    # base output dir for signal
    out_dir_base = opts.sigFile[:-len('Reco.h5')]

    for xsec, qr_run in dd.items():

        ### do bunp-hunt

        ## paths

        # distinct input dir for each cross section
        in_dir = '/eos/user/k/kiwoznia/data/QR_results/events/qr_run_'+str(dd[xsec])+'/env_run_0/poly_run_0/'
        # distinct output subdir for each cross section
        out_dir = os.path.join(out_dir_base,'xsec'+int(xsec)+'_qr'+int(qr_run))

        xsec = float(xsec)
        ysig = {}
        ypvalue = {}
        outfiles = []
        for q in quantiles:
            ysig[q] = []
            ypvalue[q] = []
            outfiles.append(open(out_dir_base+'/results_%s.txt'%q,'w'))

        ## run dijetfit
        if not opts.signi:
            cmd = "python dijetfit.py -i {in_dir} --sig {sigfile} --qcd {qcdfile} --xsec {xsec} -M {mass} --res {res} --out {out_dir}".format(in_dir=in_dir, xsec=xsec, sigfile=opts.sigFile, qcdfile=opts.qcdFile, mass=opts.mass, res=opts.sigRes, out_dir=out_dir)
            if options.correlateB == True: 
                cmd += ' -C'
            else:
                out_dir += '_noco'
            print(cmd)
            os.system(cmd)

        ## calculate significance
        for iq,q in enumerate(quantiles):

            cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{quantile}.root -m {mass} -n significance_{xsec}_{quantile}'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(mass))
            print(cmd)
            os.system(cmd)

            cmd = 'cd {out_dir} && combine -M Significance workspace_JJ_{xsec}_{quantile}.root -m {mass} -n pvalue_{xsec}_{quantile} --pvalue'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(mass))
            print(cmd)
            os.system(cmd)
                
            tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}_{quantile}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(mass)), 'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ysig[q].append(tree.limit)
            print("Xsec {}, quantile {}, significance {}".format(xsec,q,ysig[q][-1]))       
            tf.Close()

            tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}_{quantile}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec, quantile=q, mass=int(mass)), 'READ')
            tree = tf.limit
            tree.GetEntry(0)         
            ypvalue[q].append(tree.limit)        
            tf.Close()

            outfiles[iq].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=xsec,pvalue=ypvalue[q][-1],sig=ysig[q][-1]))

        # end for each quantile

        for iq, q in enumerate(quantiles): outfiles[iq].close()

    # end for each xsec

    ysig['combo'] = []
    ypvalue['combo'] = []
    outfiles.append(open(out_dir_base+'/results_final.txt','w')) # save final quantile-combined xsec scan results to signal root directory

    for xsec, qr_run in dd.items():

        out_dir = os.path.join(out_dir_base,'xsec'+int(xsec)+'_qr'+int(qr_run))
        xsec = float(xsec)

        cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n significance_{xsec}'.format(out_dir=out_dir, xsec=xsec,label='final',mass=int(mass))
        print cmd
        os.system(cmd)

        cmd = 'cd {out_dir} && combine -M Significance workspace_{xsec}_{label}.root -m {mass} -n pvalue_{xsec} --pvalue'.format(out_dir=out_dir, xsec=xsec,label='final',mass=int(mass))
        print cmd
        os.system(cmd)
            
        tf = ROOT.TFile.Open('{out_dir}/higgsCombinesignificance_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec,mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ysig['combo'].append(tree.limit)             
        print "Xsec",x,"COMBO significance",ysig['combo'][-1]        
        tf.Close()

        tf = ROOT.TFile.Open('{out_dir}/higgsCombinepvalue_{xsec}.Significance.mH{mass}.root'.format(out_dir=out_dir, xsec=xsec,mass=int(mass)),'READ')
        tree = tf.limit
        tree.GetEntry(0)             
        ypvalue['combo'].append(tree.limit)          
        tf.Close()

        outfiles[-1].write('{xsec}\t{pvalue}\t{sig}\n'.format(xsec=xsec, pvalue=ypvalue['combo'][-1], sig=ysig['combo'][-1]))  
 
    outfiles[-1].close()
    
    print ysig
    print ypvalue


    ### do GOF

