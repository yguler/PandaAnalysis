#!/bin/bash

rm -rf cache
mkdir cache

for f in ZpWW ZpTT ZpA0h QCD; do
    submit --exec convert.py --arglist partitions/${f}.txt --cache $(readlink -f cache/$f)
done
