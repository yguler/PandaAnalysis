#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20171220_common008.cfg" 
export PANDA_FLATDIR="/home/snarayan/home000/store/panda/dummy/"
mkdir -p $PANDA_FLATDIR

export SUBMIT_TMPL="skim_vbf_tmpl.py"
export SUBMIT_NAME="dummy"
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/logs/"
export SUBMIT_LOCKDIR="/data/t3serv014/snarayan/jobs/"${SUBMIT_NAME}"/locks/"
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR $SUBMIT_LOGDIR $SUBMIT_LOCKDIR


export SUBMIT_CONFIG=T2
