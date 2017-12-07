#!/bin/bash

LAST=1

while true
do
    RECENT=$(stat -c%Z $(ls -ltr $SUBMIT_LOCKDIR/* | tail -n1 | awk '{ print $9; }'))
    if [[ $RECENT > $LAST || $(expr $(date +%s) - $LAST) > 60 ]]; 
    then
        clear
        ${CMSSW_BASE}/src/PandaAnalysis/T3/bin/checkMissingFiles.py
        LAST=$(date +%s) # this takes care of the 60s refresh
        date
    fi
    sleep 5
done
