import numpy as np
import optparse
import json
import pathlib
import os


if __name__ == "__main__":

	#python main_xsec_scan.py -M 3500 --sig GtoWW35naReco.h5 --qcd qcdSigAll.h5 --res na -C

    parser = optparse.OptionParser()
    parser.add_option("-M","-M", dest="mass", type=float, default=3500., help="Injected signal mass")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcdSigAll.h5', help="QCD signalregion h5 file")
    parser.add_option("--sig","--sig", dest="sigFile", default='signal.h5', help="Signal h5 file")
    parser.add_option("--res", "--res", dest="sigRes", type="choice", choices=("na", "br"), default="na", help="resonance type: narrow [na] or broad [br]")
    parser.add_option('-C', dest="correlateB",action="store_true",help="Coorelate background shape among quantiles")
    (opts,args) = parser.parse_args()



	experiment_dict = {
		0  : 416,
		50 : 6161,
		100: 6160
	}

	# base output dir for signal
	out_dir_base = opts.sigFile[:-len('Reco.h5')]

	for i, xsec, qr_run in experiment_dict.items():

		### do dijetfit

		## paths

		# distinct input dir for each cross section
		in_dir = '/eos/user/k/kiwoznia/data/QR_results/events/qr_run_'+str(experiment_dict[xsec])+'/env_run_0/poly_run_0/'
		# distinct output subdir for each cross section
		out_dir = os.path.join(out_dir_base,'xsec'+int(xsec)+'_qr'+int(qr_run)

		## command

		cmd = "python dijetfit.py -i {in_dir} --sig {sigfile} --qcd {qcdfile} --xsec {xsec} -M {mass} --res {res} --out {out_dir}".format(in_dir=in_dir, xsec=xsec, sigfile=opts.sigFile, qcdfile=opts.qcdFile, mass=opts.mass, res=opts.sigRes, out_dir=out_dir)
        if options.correlateB == True: 
        	cmd += ' -C'
        else:
        	out_dir += '_noco'
        print(cmd)
        os.system(cmd)

        

		### do gof

