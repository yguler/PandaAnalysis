#!/usr/bin/env python

import json
import socket
from re import sub
from sys import exit
from time import clock,time
from os import system,getenv,path

import ROOT as root
from PandaCore.Tools.Misc import *
from PandaCore.Tools.Load import Load
import PandaCore.Tools.job_config as cb

sname = 'T3.job_utilities'                                     # name of this module
data_dir = getenv('CMSSW_BASE') + '/src/PandaAnalysis/data/'   # data directory
host = socket.gethostname()                                    # where we're running
IS_T3 = (host[:2] == 't3')                                     # are we on the T3?
REMOTE_READ = True                                             # should we read from hadoop or copy locally?
local_copy = bool(smart_getenv('SUBMIT_LOCALACCESS', True))    # should we always xrdcp from T2?


# global to keep track of how long things take
_stopwatch = time() 
def print_time(label):
    global _stopwatch
    now_ = time()
    PDebug(sname+'.print_time:',
           '%.1f s elapsed performing "%s"'%((now_-_stopwatch),label))
    _stopwatch = now_


# convert an input name to an output name
def input_to_output(name):
    if 'input' in name:
        return name.replace('input','output')
    else:
        return 'output_' + name.split('/')[-1]


# find data and bring it into the job somehow 
#  - if local_copy and it exists locally, then:
#    - if REMOTE_READ, read from hadoop
#    - else copy it locally
#  - else xrdcp it locally
def copy_local(long_name):
    full_path = long_name
    PInfo(sname,full_path)

    panda_id = long_name.split('/')[-1].split('_')[-1].replace('.root','')
    input_name = 'input_%s.root'%panda_id
    copied = False

    # if the file is cached locally, why not use it?
    local_path = full_path.replace('root://xrootd.cmsaf.mit.edu/','/mnt/hadoop/cms')
    PInfo(sname+'.copy_local','Local access is configured to be %s'%('on' if local_copy else 'off'))
    if local_copy and path.isfile(local_path): 
        # apparently SmartCached files can be corrupted...
        ftest = root.TFile(local_path)
        if ftest and not(ftest.IsZombie()):
            PInfo(sname+'.copy_local','Opting to read locally')
            if REMOTE_READ:
                return local_path
            else:
                cmd = 'cp %s %s'%(local_path, input_name)
                PInfo(sname+'.copy_local',cmd)
                system(cmd)
                copied = True

    if not copied:
        cmd = "xrdcp %s %s"%(full_path,input_name)
        PInfo(sname+'.copy_local',cmd)
        system(cmd)
        copied = True
            
    if path.isfile(input_name):
        PInfo(sname+'.copy_local','Successfully copied to %s'%(input_name))
        return input_name
    else:
        PError(sname+'.copy_local','Failed to copy %s'%input_name)
        return None


# wrapper around rm -f. be careful!
def cleanup(fname):
    ret = system('rm -f %s'%(fname))
    if ret:
        PError(sname+'.cleanup','Removal of %s exited with code %i'%(fname,ret))
    else:
        PInfo(sname+'.cleanup','Removed '+fname)
    return ret


# wrapper around hadd
def hadd(good_inputs, output='output.root'):
    good_outputs = ' '.join([input_to_output(x) for x in good_inputs])
    cmd = 'hadd -f ' + output + ' ' + good_outputs
    print cmd
    ret = system(cmd)    
    if not ret:
        PInfo(sname+'.hadd','Merging exited with code %i'%ret)
    else:
        PError(sname+'.hadd','Merging exited with code %i'%ret)


# remove any irrelevant branches from the final tree.
# this MUST be the last step before stageout or you 
# run the risk of breaking something
def drop_branches(to_drop=None, to_keep=None):
    if not to_drop and not to_keep:
        return 0

    if to_drop and to_keep:
        PError(sname+'.drop_branches','Can only provide to_drop OR to_keep')
        return 0

    f = root.TFile('output.root','UPDATE')
    t = f.FindObjectAny('events')
    n_entries = t.GetEntriesFast() # to check the file wasn't corrupted
    if to_drop:
        if type(to_drop)==str:
            t.SetBranchStatus(to_drop,False)
        else:
            for b in to_drop:
                t.SetBranchStatus(b,False)
    elif to_keep:
        t.SetBranchStatus('*',False)
        if type(to_keep)==str:
            t.SetBranchStatus(to_keep,True)
        else:
            for b in to_keep:
                t.SetBranchStatus(b,True)
    t_clone = t.CloneTree()
    f.WriteTObject(t_clone,'events','overwrite')
    f.Close()

    # check that the write went okay
    f = root.TFile('output.root')
    if f.IsZombie():
        PError(sname+'.drop_branches','Corrupted file trying to drop '+str(to_drop))
        return 1 
    t_clone = f.FindObjectAny('events')
    if (n_entries==t_clone.GetEntriesFast()):
        return 0
    else:
        PError(sname+'.drop_branches','Corrupted tree trying to drop '+str(to_drop))
        return 2


# stageout a file (e.g. output or lock)
#  - if IS_T3, execute a simple mv
#  - else, use lcg-cp
# then, check if the file exists:
#  - if IS_T3, use os.path.isfile
#  - else, use lcg-ls
def stageout(outdir,outfilename,infilename='output.root'):
    if IS_T3:
        mvargs = 'mv $PWD/%s %s/%s'%(infilename,outdir,outfilename)
    else:
        mvargs = 'lcg-cp -v -D srmv2 -b file://$PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(infilename,outdir,outfilename)
        #mvargs = 'gfal-copy $PWD/%s srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(infilename,outdir,outfilename)
    PInfo(sname+'.stageout',mvargs)
    ret = system(mvargs)
    if not ret:
        PInfo(sname+'.stageout','Move exited with code %i'%ret)
    else:
        PError(sname+'.stageout','Move exited with code %i'%ret)
        return ret
    if IS_T3:
        if not path.isfile('%s/%s'%(outdir,outfilename)):
            PError(sname+'.stageout','Output file is missing!')
            ret = 1
    else:
        lsargs = 'lcg-ls -v -D srmv2 -b srm://t3serv006.mit.edu:8443/srm/v2/server?SFN=%s/%s'%(outdir,outfilename)
        PInfo(sname+'.stageout',lsargs)
        ret = system(lsargs)
        if ret:
            PError(sname+'.stageout','Output file is missing!')
    return ret


# write a lock file, based on what succeeded,
# and then stage it out to a lock directory
def write_lock(outdir,outfilename,processed):
    outfilename = outfilename.replace('.root','.lock')
    flock = open(outfilename,'w')
    for k,v in processed.iteritems():
        flock.write(v+'\n')
    flock.close()
    stageout(outdir,outfilename,outfilename)


# classify a sample based on its name
def classify_sample(full_path, isData):
    if not isData:
        if any([x in full_path for x in ['Vector_','Scalar_']]):
            return root.kSignal
        elif any([x in full_path for x in ['ST_','ZprimeToTT']]):
            return root.kTop
        elif 'EWKZ2Jets' in full_path:
            return root.kZEWK
        elif 'EWKW' in full_path:
            return root.kWEWK
        elif 'ZJets' in full_path or 'DY' in full_path:
            return root.kZ
        elif 'WJets' in full_path:
            return root.kW
        elif 'GJets' in full_path:
            return root.kA
        elif 'TTJets' in full_path or 'TT_' in full_path:
            return root.kTT
    return root.kNoProcess


# read a CERT json and add it to the skimmer
def add_json(skimmer, json_path):
    with open(json_path) as jsonFile:
        payload = json.load(jsonFile)
        for run_str,lumis in payload.iteritems():
            run = int(run_str)
            for l in lumis:
                skimmer.AddGoodLumiRange(run,l[0],l[1])


# some common stuff that doesn't need to be configured
def run_PandaAnalyzer(skimmer, isData, input_name):
    # read the inputs
    try:
        fin = root.TFile.Open(input_name)
        tree = fin.FindObjectAny("events")
        weight_table = fin.FindObjectAny('weights')
        hweights = fin.FindObjectAny("hSumW")
    except:
        PError(sname+'.run_PandaAnalyzer','Could not read %s'%input_name)
        return False # file open error => xrootd?
    if not tree:
        PError(sname+'.run_PandaAnalyzer','Could not recover tree in %s'%input_name)
        return False
    if not hweights:
        PError(sname+'.run_PandaAnalyzer','Could not recover hweights in %s'%input_name)
        return False
    if not weight_table:
        weight_table = None

    output_name = input_to_output(input_name)
    skimmer.SetDataDir(data_dir)
    if isData:
        add_json(skimmer, data_dir+'/certs/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt')

    rinit = skimmer.Init(tree,hweights,weight_table)
    if rinit:
        PError(sname+'.run_PandaAnalyzer','Failed to initialize %s!'%(input_name))
        return False 
    skimmer.SetOutputFile(output_name)

    # run and save output
    skimmer.Run()
    skimmer.Terminate()

    ret = path.isfile(output_name)
    if ret:
        PInfo(sname+'.run_PandaAnalyzer','Successfully created %s'%(output_name))
        return output_name 
    else:
        PError(sname+'.run_PandaAnalyzer','Failed in creating %s!'%(output_name))
        return False


# main function to run a skimmer, customizable info 
# can be put in fn
def main(to_run, processed, fn):
    print_time('loading')
    for f in to_run.files:
        input_name = copy_local(f)
        print_time('copy %s'%input_name)
        if input_name:
            success = fn(input_name,(to_run.dtype!='MC'),f)
            print_time('analyze %s'%input_name)
            if success:
                processed[input_name] = f
            if input_name[:5] == 'input': # if this is a local copy
                cleanup(input_name)
            print_time('remove %s'%input_name)
    
    if len(processed)==0:
        PWarning(sname+'.main', 'No successful outputs!')
        exit(1)


