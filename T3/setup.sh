#!/bin/bash

export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
export PANDA_CFG="http://t3serv001.mit.edu/~snarayan/histcatalog/20171109_tag.cfg" 
export PANDA_FLATDIR="${HOME}/home000/store/panda/v_005_deep0/"
mkdir -p $PANDA_FLATDIR

export SUBMIT_TMPL="skim_deep_tmpl.py"
export SUBMIT_NAME="v_005_deep0"
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/logs/"
export SUBMIT_OUTDIR="/data/t3serv014/snarayan/panda/"${SUBMIT_NAME}"/batch/"
mkdir -p $SUBMIT_WORKDIR $SUBMIT_OUTDIR/locks/ $SUBMIT_LOGDIR

export SUBMIT_LOCALACCESS=1
