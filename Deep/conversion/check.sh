#!/bin/bash

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in Top QCD Higgs W; do
    echo $f
    check --cache cache/$f $@
done
