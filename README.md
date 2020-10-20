# Dijet fit with shapes for anomaly detection with VAE

### Prerequisites

Higgs combine tools is needed, either standalone or from cmssw. See [here](http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/) for instructions to get latest version.

### 1-category bump hunt

Run the 1-category and N-category s+b fit with injected signal of a chosen cross-section `--xsec`. The signal mass hypothesis can also be configured with `-M`. Signal and background
h5 files are taken from the input folder set with `-i`. The `--sig` and `--qcd` should be relative to that input folder.

```
python dijetfit.py --xsec {XSEC} -M {MASS} -i {INPUTDIR} --sig {SIGNALPATH} --qcd {QCDPATH}
```

### Automatic scan of p-value/significance vs cross section

```
python run_dijetfit.py --run --i {INPUTDIR} -M {MASS} -i {INPUTDIR} --sig {SIGNALPATH} --qcd {QCDPATH}
```

This first makes the workspaces input to combine and then computes expected significance with pre-fit Asimov dataset per each signal xsec hypothesis.
The cross section range can be change in the `files_count.json` file. The results are saved in txt files for each quantile and for the combination. If you have run the full scan and you just
want to plot from the txt files then remove option `--run`.
