#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/test_t2.cfg" 
export PANDA_FLATDIR="/home/snarayan/home000/store/panda/v_008_v0/"
mkdir -p $PANDA_FLATDIR

export SUBMIT_TMPL="skim_gghbb_tmpl.py"
export SUBMIT_NAME="v_008_v0"
export SUBMIT_WORKDIR="/data/t3home000/snarayan/panda/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3home000/snarayan/panda/panda/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR/locks/ $SUBMIT_LOGDIR


# testing
export SUBMIT_CONFIG=T2
