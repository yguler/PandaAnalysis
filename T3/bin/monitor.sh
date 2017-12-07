#!/bin/bash

LASTMODIFIED=1

while true
do
    RECENTLYMODIFIED=$(stat -c%Z $(ls -ltr $SUBMIT_LOCKDIR/* | tail -n1 | awk '{ print $9; }'))
    if [[ $RECENTLYMODIFIED > $LASTMODIFIED ]]; 
    then
        clear
        ${CMSSW_BASE}/src/PandaAnalysis/T3/bin/checkMissingFiles.py
        LASTMODIFIED=$RECENTLYMODIFIED
    fi
    sleep 5
done
