#!/bin/bash

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in cache/*; do
    echo $f
    check --cache $f $@
done
