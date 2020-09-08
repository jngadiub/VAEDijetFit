# Dijet fit with shapes for anomaly detection with VAE

### 1-category bump hunt

Run the 1-category s+b fit with injected signal of a chosen cross-section `--xsec`.Chose the region to fit with `--index` where 0=inclusive,
1=accepted 2=rejected. The signal mass hypthesis can also be configured with `-M`. Name and path of signa/background h5 files can be set with `--qcd` and `--sig`. 

```
python dijetfit.py --index {REGION} --xsec {XSEC} -M {MASS}
```

### 2-category bump hunt

The 2-category s+b fit is run in two steps:


1. First make the control region datacard with the rejected sample with the additional `--twoCatFit` option:

```
python dijetfit.py --index 2 --xsec {XSEC} -M {MASS} --twoCatFit
```

2. Then make the signal region datacard which is combined with the control region datacard produced in step 1

```
python dijetfit.py --index 1 --xsec {XSEC} -M {MASS} --twoCatFit
```

### Automatic scan of p-value/significance vs cross section

```
python run_dijetfit.py 1
```

This runs the scan for the 1-category bump hunt for both the inclusive and accepted samples. After that it runs the 2-category bump hunt. The cross
section range and mass can be changed in the script. The results are saved in txt files (they should be three). If you have run the full scan and you just
want to plot from the txt files then set first argument to 0:

```
python run_dijetfit.py 0
```

