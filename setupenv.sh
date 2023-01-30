#cd /afs/cern.ch/user/k/kiwoznia/CMSSW_10_2_3/src
#cd /afs/cern.ch/user/k/kiwoznia/CMSSW_11_3_4/src
cd /afs/cern.ch/user/k/kiwoznia/CMSSW_10_2_13/src # DO NOT TOUCH OTHERWISE SEGFAULT ARMAGEDDON
source /cvmfs/cms.cern.ch/cmsset_default.sh
eval `scramv1 runtime -sh`
#cmsenv
export PYTHONPATH=
export PYTHON27PATH=$PYTHON27PATH:/afs/cern.ch/user/k/kiwoznia/.local/lib/python2.7/site-packages
cd -

########################################################################################

# for the GOF test run 
# source /cvmfs/sft.cern.ch/lcg/views/LCG_98python3/x86_64-centos7-gcc10-opt/setup.sh
#