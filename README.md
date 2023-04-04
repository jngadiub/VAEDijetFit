# Dijet fit with shapes for anomaly detection with VAE

### Prerequisites

Higgs combine tools is needed, either standalone or from cmssw. See [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) for instructions to get latest version.

### N-category bump hunt for a signal cross section hypothesis

Run the individual 1-category and N-category s+b fit with manually injected signal of a chosen cross-section `--xsec`. The signal mass hypothesis must also be configured with `-M`. Signal and background h5 files are taken from the input folder set with `-i`. The `--sig` and `--qcd` should be relative to that input folder. 

For implementation of full correlation of background shapes and normalizations among quantiles add option `-C`. If you want to run combine FitDiagnostics and Significance to get observed signal strenght and significance add option `-R`.

Example command used for [latest CASE ANv7](http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2020_051_v7.pdf):

```
python dijetfit.py -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ --sig signal_XToYYprimeTo4Q_MX2000_MY80_MYprime170_narrowReco.h5 --qcd bkg_XToYYprimeTo4Q_MX2000_MY80_MYprime170_narrowReco_0.25.h5 --xsec 0.25 -M 2000.0 --out test -C -R
```

Several quantiles configuration can be tested by setting `--config`. The default value is now `4` corresponding to the latest chosen configuration with only 3 quantiles (90%,95%,99%).

### Automatic scan of p-value/significance vs cross section

The options are similar to the main `dijet.py` script with a few additional ones. With option `--init` the expected significance from post-fit Asimov is computed. Remove the option for manual signal injection and oberved significance computation.

The cross section range can be changed in the `files_count.json` file. The results are saved in txt files for each quantile and for the combination. If you have run the full scan and you just want to plot from the txt files then remove option `--run`.

Example command for results in [latest CASE ANv7](http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2020_051_v7.pdf):

```
python run_dijetfit.py --i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 --qcd bkg.h5 -M 2000 -n 1 --run -C
```

### Bias test

The bias test is computed by first making a workspace with background-only data but with meaningful signal normalization running the main `dijetfit.py` script and then by running the post-Asimov expected signicance with combine for many toys over several condor jobs. You can set the number of condor jobs and the number of toys per job with options `-j` and `-t` respectively, Once the jobs are finished remove option `-s` for plotting the results.

Example command for results in [latest CASE ANv7](http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2020_051_v7.pdf):

```
python submit-jobs.py -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec0_init --qcd bkg.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 0.0 --exp_sig 0 -s

python submit-jobs.py -j 400 -t 10 -i /afs/cern.ch/work/i/izoi/public/forJennifer/run_41025/ -w bias_test_xsec1_init --qcd bkg.h5 --sig XToYYprimeTo4Q_MX2000_MY80_MYprime170 -M 2000 --xsec 1.0 --exp_sig 0.25 -s
```

Note: use `--xsec 0.0` if you run the bias test for zero expected signal (i.e. `--exp_sig 0`). Otherwise use `--xsec 1.0` for bias test with non-zero expected signal. This is done so that the number used for the `--expectSignal` is equal to the cross section.