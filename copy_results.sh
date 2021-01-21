#!/bin/bash
echo "creating directory for quantile $1 and strategy $2"
export RES_DIR="$1_$2_results"
mkdir ${RES_DIR}
mv *.root ${RES_DIR}
mv *.txt ${RES_DIR}
mv pvalue.png ${RES_DIR}

