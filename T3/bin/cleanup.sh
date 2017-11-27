#!/bin/bash

cd $SUBMIT_OUTDIR
pwd

sleep 2
rm -f *root *npz locks/* &

cd -
