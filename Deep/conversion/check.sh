#!/bin/bash

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in cache/*; do
    PInfo -n check.sh $f
    check --cache $f $@
done
