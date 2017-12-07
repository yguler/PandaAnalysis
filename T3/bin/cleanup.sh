#!/bin/bash

echo 'Cleaning up!'
echo $SUBMIT_OUTDIR
echo $SUBMIT_LOCKDIR

sleep 2

rm -rf $SUBMIT_OUTDIR/* &
rm -rf $SUBMIT_LOCKDIR/*
