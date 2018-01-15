# Analysis on T3/T2 Condor

Throughout this readme, I will refer to several user-defined environment variables. 
These are typically defined in `T3/setup.sh`, but the user can define things elswhere.

## Cataloging inputs

```bash
./bin/catalogT2Prod.py --outfile /path/to/config.cfg [ --include datasets to include ] [ --exclude datasets to skip ] [ --force ] [--smartcache]
```

I recommend you put the config in a web-accessible place for use in later steps. For example:
```bash
./catalogT2Prod.py --force --outfile ~/www/eoscatalog/$(date +%Y%m%d).cfg --include TT --exclude TTbarDM --smartcache
```

The above command will do the following things:

- It will only check datasets that contain `TT` and do not contain `TTbarDM` in the dataset's nickname

- If a dataset is not found in `PandaCore.Tools.process`, it will guess the nickname of the dataset, give it a xsec of 1, and write it to the catalog (`--force`)

- If the file does not exist locally on the T3, a smartcache request will be made

- The output will be a timestamped web-facing config file


## Building the work environment

First, make sure the following environment variables are defined (some examples shown):
```bash
export PANDA_CFG="http://snarayan.web.cern.ch/snarayan/eoscatalog/20170127.cfg"  # location of config file from previous section
export SUBMIT_TMPL="skim_merge_tmpl.py"  # name of template script in T3/inputs
export SUBMIT_NAME="v_8024_2_0"  # name for this job
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/work/"  # staging area for submission
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/logs/"  # log directory
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"  # location of unmerged files
export PANDA_FLATDIR="${HOME}/home000/store/panda/v_8024_2_0/"   # merged output
```

`T3/inputs/$SUBMIT_TMPL` should be the skimming configuration you wish to run your files through. 
This will be replaced by a `sed` with the output directory path before submission.

The inputs are then built by doing:
```bash
./buildMergedInputs.sh -t [-n number of files per job]
```
The default if `-n` is not provided is 20. 
I recommend you limit the number of jobs to less than 700, which means typically a value of 40-50.
Submitting more jobs won't make them run any faster.

## Submitting and re-submitting

To submit jobs, simply do

```bash 
python submit.py
```
To check the status of your jobs, simply do:
```bash
./checkMissingFiles.py  [--silent] [--force] [--nfiles NFILES]
```
Note that the above command overwrites `$SUBMIT_WORKDIR/local.cfg`, with the intention of preparing it for resubmission.
The file will be recreated as a configuration to rerun files that are not present in the output and not running.
The option `--silent` will skip the per-sample breakdown.
The option `--nfiles` will repackage `local.cfg` into a different number of files per job.
The option `--force` will re-catalog files that are incomplete, not just missing.

To resubmit missing files, simply do
```bash
python submit.py
```
In the case that you are using the `--force` option, make sure you have no running jobs before resubmitting, or you may end up with duplicated outputs.


## Merging

Make sure `$PANDA_FLATDIR` exists. Then, go into `T3/merging` and do:
```bash
./merge.py TTbar_Powheg
```
to merge the Powheg TT sample, for example. 
To merge en-masse (e.g. many many signal outputs), you can do something like:
```bash
submit --exec merge.py --arglist list_of_signals.txt
```
This assumes that you have `PandaCore/bin` in your `$PATH`
