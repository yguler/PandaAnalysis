#!/bin/bash                                                                                                                                                                                                 
export PATH=${PATH}:${CMSSW_BASE}/src/PandaCore/bin/

#submission number
export SUBMIT_NAME="v_8026_monoj_v2"
#scratch space
export scratch_area="/uscms_data/d3"
export PANDA="${CMSSW_BASE}/src/PandaAnalysis"
#cfg file
export PANDA_CFG="http://t3serv001.mit.edu/~mcremone/eoscatalog/20170905.cfg"
#skim
export SUBMIT_TMPL="skim_monoj_tmpl.py"
#panda's 
export PANDA_FLATDIR="${scratch_area}/${USER}/panda/"${SUBMIT_NAME}"/flat/"
export SUBMIT_OUTDIR="/store/user/${USER}/panda/"${SUBMIT_NAME}"/batch/"
#condor's
export SUBMIT_WORKDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/work/"
export SUBMIT_LOGDIR="${scratch_area}/${USER}/condor/"${SUBMIT_NAME}"/logs/"
mkdir -p $PANDA_FLATDIR $SUBMIT_WORKDIR $SUBMIT_LOGDIR
eosmkdir -p $SUBMIT_OUTDIR

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo ""
echo ""
cat << "EOF"
  _____        _   _ _____                      ______ _   _          _          _____  __  __ 
 |  __ \ /\   | \ | |  __ \   /\        ____   |  ____| \ | |   /\   | |        |  __ \|  \/  |
 | |__) /  \  |  \| | |  | | /  \      / __ \  | |__  |  \| |  /  \  | |  ______| |  | | \  / |
 |  ___/ /\ \ | . ` | |  | |/ /\ \    / / _` | |  __| | . ` | / /\ \ | | |______| |  | | |\/| |
 | |  / ____ \| |\  | |__| / ____ \  | | (_| | | |    | |\  |/ ____ \| |____    | |__| | |  | |
 |_| /_/    \_\_| \_|_____/_/    \_\  \ \__,_| |_|    |_| \_/_/    \_\______|   |_____/|_|  |_|
                                       \____/                                                  
EOF
echo ""
echo "Checking ENV path"
echo "======================================================================="

for path in $PANDA_FLATDIR /eos/uscms${SUBMIT_OUTDIR} $SUBMIT_WORKDIR $SUBMIT_LOGDIR
do
if [ -e $path ];then
echo "Path : ${path} is properly set"
else
echo "Path : ${path} does not exist, please fix it."
fi
done
echo "======================================================================"
echo "INFO"
echo "======================================================================"
echo "Submit Name  = ${SUBMIT_NAME}"
echo "cfg selected = ${PANDA_CFG}"
echo "submit tmpl  = ${SUBMIT_TMPL}"
echo "======================================================================"
echo ""
echo ""

#private production                                                                                                                                                       
#export PRIVATE_LOGDIR="${HOME}/cms/logs/monotop_private_panda/"
#export PRIVATE_PRODDIR="${HOME}/cms/hist/monotop_private_pandatree/"
#export PRIVATE_CFGDIR="${HOME}/cms/condor/monotop_private_panda/"

# fitting                                                                                                                                                                                                   
#export PANDA_FIT=/data/t3serv014/snarayan/CMSSW_7_4_7/
#export PANDA_XSECS=/home/snarayan/cms/cmssw/analysis/MonoTop_Xsec/
#export PANDA_XSECS=/eos/uscms/store/user/shoh/zprime_cross_section/
#export PANDA_PROD=/eos/uscms/store/user/shoh/miniaod/
#export PANDA_FITTING=${PANDA_FLATDIR}/fitting/
#mkdir -p $PANDA_FITTING/scans/ $PANDA_FITTING/logs/


