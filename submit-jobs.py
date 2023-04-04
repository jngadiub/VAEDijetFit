import os, sys
import numpy as np
import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import ROOT
import optparse
from scipy.stats import norm

def makeSubmitFileCondor(exe,jobname,jobflavour):
    print "make options file for condor job submission "
    submitfile = open("submit.sub","w")        
    submitfile.write("should_transfer_files = YES\n")
    submitfile.write("when_to_transfer_output = ON_EXIT\n")
    submitfile.write('transfer_output_files = ""\n')
    submitfile.write("executable  = "+exe+"\n")
    submitfile.write("arguments             = $(ClusterID) $(ProcId)\n")    
    submitfile.write("output                = "+jobname+".$(ClusterId).$(ProcId).out\n")
    submitfile.write("error                 = "+jobname+".$(ClusterId).$(ProcId).err\n")
    submitfile.write("log                   = "+jobname+".$(ClusterId).log\n")
    submitfile.write('+JobFlavour           = "'+jobflavour+'"\n')
    submitfile.write("queue")
    submitfile.close() 

if __name__ == "__main__": 

    #python submit-jobs.py -s -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec0_init --qcd bkg.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 0.0 --exp_sig 0
    #python submit-jobs.py -s -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec1_init --qcd bkg.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 1.0 --exp_sig 0.25
    #python submit-jobs.py -s -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec1_init --qcd bkg.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 1.0 --exp_sig 0.5
    #python submit-jobs.py -s -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec1_init --qcd bkg_XToYYprimeTo4Q_MX2000_MY80_MYprime170_narrowReco_0.25.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 1.0 --exp_sig 0.25

    parser = optparse.OptionParser()
    parser.add_option("-s", "--submit_jobs", dest="submit_jobs", action="store_true", default=False, help="Submit jobs")
    parser.add_option("-j","--jobs",dest="jobs",type=int,default=2,help="Number of jobs")
    parser.add_option("-t","--toys",dest="toys",type=int,default=2,help="Number of toys per job")
    parser.add_option("--sig","--sig",dest="signal",type=str,default='XToYYprimeTo4Q_MX2000_MY80_MYprime170',help="Signal name")
    parser.add_option("--add_text","--add_text",dest="add_text",type=str,default='No signal injected',help="Additional text for plot")
    parser.add_option("--exp_sig","--exp_sig",dest="exp_sig",type=float, default=0,help="Expected signal for combine")
    parser.add_option("--qcd","--qcd", dest="qcdFile", default='qcd.h5', help="QCD h5 file")
    parser.add_option("--xsec","--xsec",dest="xsec",type=float,default=0.0,help="Injected signal cross section in pb")
    parser.add_option("-M","-M",dest="mass",type=float,default=3500.,help="Injected signal mass")
    parser.add_option("-C","--config",dest="config",type=str,default='-C --config 4',help="Analysis configuration (see dijetfit.py)")
    parser.add_option("-i","--inputDir",dest="inputDir",type=str,default='./',help="Directory with all quantiles h5 files")
    parser.add_option("-w","--workspace_dir",dest="workspace_dir",type=str,default='./',help="Directory with the initial workspaces (local path)")

    (options,args) = parser.parse_args()
    print(options)

    submit = options.submit_jobs
    jobs = options.jobs
    toys = options.toys
    indir = options.inputDir
    workspace_dir = os.path.join(os.getcwd(),options.workspace_dir)
    signal=options.signal
    mass = options.mass
    xsec = options.xsec
    exp_signal = options.exp_sig
    qcdfile = options.qcdFile
    add_text = options.add_text
    add_text = add_text.replace('\sigma','$\sigma$')

    jdir = "bias_test_injQR_{SIG}_xsec{EXPSIG}".format(SIG=signal,EXPSIG=exp_signal)
    outdir = os.getcwd()
    quantiles = ['_q90', '_q95', '_q99', '_total', '_final']

    if submit==True:

      cmd = "python dijetfit.py --xsec {XSEC} -M {MASS} -i {INDIR} --sig signal_{SIGNAL}_narrowReco.h5 --qcd {QCDFILE} --out {OUTDIR} -C --config 4 --run_toys".format(INDIR=indir,OUTDIR=workspace_dir,MASS=mass,SIGNAL=signal,XSEC=xsec,QCDFILE=qcdfile)
      print(cmd)
      os.system(cmd)

      for j in range(0,jobs):

        jobdir = jdir+"_job{ID}".format(ID=j+1)
        os.mkdir(jobdir)  
        os.chdir(jobdir)
        print(jobdir)

        outdir=os.path.join(os.getcwd())       
        my_cmd = ""
        for q in quantiles:
          my_cmd += "combine -M FitDiagnostics -m {MASS} -n {Q} --toysFreq -t {TOYS} --rMin -100 -s -1 --expectSignal={EXPSIG} {INDIR}/workspace_JJ_{XSEC}{Q}.root --out {OUTDIR}\n".format(TOYS=toys,INDIR=workspace_dir,Q=q,OUTDIR=outdir,MASS=mass,EXPSIG=exp_signal,XSEC=xsec*1000.)

        with open('job.sh', 'w') as fout:
          fout.write("#!/bin/sh\n")
          fout.write("echo\n")
          fout.write("echo\n")
          fout.write("echo 'START---------------'\n")
          fout.write("echo 'WORKDIR ' ${PWD}\n")
          fout.write("source /afs/cern.ch/cms/cmsset_default.sh\n")
          fout.write("cd "+str(outdir)+"\n")
          fout.write("cmsenv\n")
          fout.write("export X509_USER_PROXY=$1\n")
          fout.write("echo $X509_USER_PROXY\n")
          fout.write("voms-proxy-info -all\n")
          fout.write("voms-proxy-info -all -file $1\n")
          fout.write("%s\n"%(my_cmd)) 
          fout.write("echo 'STOP---------------'\n")
          fout.write("echo\n")
          fout.write("echo\n")
    
        os.system("chmod 755 job.sh")      

        ###### sends job ######
        makeSubmitFileCondor("job.sh","job","microcentury")
        os.system("condor_submit submit.sub")
        print "job nr " + str(j+1) + " submitted"

        os.chdir("../")

      sys.exit()

    #quantiles = ['_q0', '_q30', '_q50', '_q70', '_q90', '_total', '']
    quantiles = ['_q90', '_q95', '_q99', '_total', '']
    results = np.zeros((len(quantiles),jobs*toys))
    errLow = np.zeros((len(quantiles),jobs*toys))
    errHigh = np.zeros((len(quantiles),jobs*toys))
    err = np.zeros((len(quantiles),jobs*toys))

    for i,q in enumerate(quantiles):

      print("Quantile",q)

      ind = 0
      for j in range(jobs):

        if j%100==0: print(j)

        jobdir = jdir+"_job{ID}".format(ID=j+1)
        if not os.path.exists(jobdir): continue

        for f in os.listdir(jobdir):
          infile = ''
          if q in f and 'fitDiagnostics' in f:
            infile = os.path.join(jobdir,f)
            break
        if os.path.exists(infile):
          try:
            tf = ROOT.TFile(infile,'READ')
            t = tf.tree_fit_sb            
          except:
            continue   
          for toy in range(toys):    
            t.GetEntry(toy)
            results[i][ind]=t.r
            err[i][ind]=t.rErr
            errLow[i][ind]=t.rLoErr
            errHigh[i][ind]=t.rHiErr
            ind+=1
          tf.Close()
        else: print "MISSING FILE FOR JOB",j+1,infile

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    for i,q in enumerate(quantiles):

      print(q)
      
      r = results[i,results[i,:] != 0]
      rErr = err[i,results[i,:] != 0]
      rErrHigh = errHigh[i,results[i,:] != 0]
      rErrLow = errLow[i,results[i,:] != 0]
      rErr[r<0] = rErrHigh[r<0]
      rErr[r>0] = rErrLow[r>0]
      r = r[rErr!=0]
      rErr = rErr[rErr!=0]
      pull = (r-exp_signal)/rErr
      print(r.shape,rErr.shape,pull.shape)   

      fig,ax = plt.subplots(figsize=(6, 6))

      mean,std=norm.fit(pull)

      font = {'size':12, 'weight': 'bold'}
      textstr = '\n'.join((
      signal,
      add_text))
      ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontdict=font, verticalalignment='top')

      textstr = '\n'.join((
      "Gaussian fit:",  
      'mean = %.2f' % (mean, ),
      'std = %.2f' % (std, )))
      ax.text(0.05, 0.75, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top')

      bins_ = np.arange(-5, 5, 0.2)
      ax.hist(pull, bins=bins_, normed=True, histtype='step')
      ax.plot(bins_, norm.pdf(bins_, mean, std), linewidth=2)
      ax.set_ylim(0,ax.get_ylim()[1]*1.20)
      ax.set_xlabel("$(r_{fit}-r_{exp})/\sigma_{r}$",fontsize=20)
      fig.tight_layout()
      fig.savefig('figures/bias%s.png'%q)    


