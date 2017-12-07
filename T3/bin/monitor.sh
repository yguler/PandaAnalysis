#!/bin/bash

LAST=1

while true
do
    RECENT=$(stat -c%Z $(ls -ltr $SUBMIT_LOCKDIR/* | tail -n1 | awk '{ print $9; }'))
    if [[ $RECENT > $LAST ]]; 
    then
        clear
        ${CMSSW_BASE}/src/PandaAnalysis/T3/bin/checkMissingFiles.py
        LAST=$RECENT
    fi
    sleep 5
done
