#!/bin/bash 

echo MET SingleElectron TTbar ZtoNuNu ZJets WJets SingleTop QCD Diboson | xargs -n 1 -P 10 python merge.py
#echo Diboson | xargs -n 1 -P 10 python merge.py
#echo MET SinglePhoton SingleElectron TTbar ZJets GJets WJets SingleTop QCD Diboson TTbar_isrup TTbar_isrdown TTbar_FXFX TTbar_tuneup TTbar_tunedown | xargs -n 1 -P 10 python merge.py
#echo MET SingleElectron SinglePhoton | xargs -n 1 -P 5 python merge.py
while read p; do
  echo $p  | xargs -n 1 -P 10 python merge.py
done <monoh.txt

