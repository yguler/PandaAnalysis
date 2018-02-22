#!/bin/bash

executable=$1

if [[ "$executable" == "" ]]; then 
    echo "Usage: $0 <executable>" >&2
    exit 1
fi

rm -rf cache
mkdir cache

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in Top_lo QCD; do
#for f in Top_lo; do
    submit --exec $executable --arglist partitions/${f}.txt --cache $(readlink -f cache/$f) 
done
