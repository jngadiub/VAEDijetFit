#!/bin/bash
echo "running dijet fit for quantile $1 and strategy $2"
export BG_FILE="/eos/project/d/dshep/TOPCLASS/DijetAnomaly/VAE_results/run_101/sample_results/$2/$1/qcd_sqrtshatTeV_13TeV_PU40_ALL_parts/qcd_sqrtshatTeV_13TeV_PU40_ALL_reco.h5"
export SIG_FILE="/eos/project/d/dshep/TOPCLASS/DijetAnomaly/VAE_results/run_101/sample_results/$2/$1/RSGraviton_WW_NARROW_13TeV_PU40_3.5TeV_parts/RSGraviton_WW_NARROW_13TeV_PU40_3.5TeV_reco.h5"
python run_dijetfit.py 1
