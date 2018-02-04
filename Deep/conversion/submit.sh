#!/bin/bash

executable=$1

rm -rf cache
mkdir cache

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in Top_lo QCD Higgs; do
    submit --exec $executable --arglist partitions/${f}.txt --cache $(readlink -f cache/$f) 
done
