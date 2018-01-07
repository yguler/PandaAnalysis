#!/usr/bin/env python

from os import getenv,path,popen,system
from PandaCore.Tools.job_management import *
import subprocess
import sys
import argparse
from glob import glob
from re import sub
import cPickle as pickle
from itertools import chain 
from query import query
from time import time, sleep, strftime
import curses
from submit import main as resubmit

lockdir = getenv('SUBMIT_LOCKDIR')
workdir = getenv('SUBMIT_WORKDIR')
parser = argparse.ArgumentParser(description='check missing files')
parser.add_argument('--infile',type=str,default=None)
parser.add_argument('--outfile',type=str,default=None)
parser.add_argument('--lockdir',type=str,default=lockdir)
parser.add_argument('--force',action='store_true')
parser.add_argument('--nfiles',type=int,default=-1)
parser.add_argument('--width',type=int,default=None)
parser.add_argument('--silent',action='store_true')
parser.add_argument('--verbose',action='store_true')
parser.add_argument('--monitor',type=int,default=None)
parser.add_argument('--resubmit',action='store_true')
args = parser.parse_args()
lockdir = args.lockdir

if args.monitor and args.verbose:
    args.verbose = False

if not args.infile:
    args.infile = workdir+'/local_all.cfg'
if not args.outfile:
    args.outfile = workdir+'/local.cfg'

columns = int(popen('stty size', 'r').read().split()[-1])
rows = int(popen('stty size', 'r').read().split()[0])
if not args.width:
    WIDTH = (columns-90)/2
else:
    WIDTH = args.width
header = ('%%-%is'%(WIDTH))%('Sample')
header += ('%%-%is'%(WIDTH+2))%('Progress')
header += ' %10s %10s %10s %10s %10s %10s'%('Total','T3', 'T2','Idle','Missing','Done')

colors = {
        'green' : 42,
        'blue' : 44,
        'cyan' : 46,
        'grey' : 47, 
        'red' : 41,
        }

last_lock = 1
last_check = 1

def init_colors():
    curses.start_color()
    curses.init_pair(colors['green'], curses.COLOR_WHITE, curses.COLOR_GREEN)
    curses.init_pair(colors['blue'], curses.COLOR_WHITE, curses.COLOR_BLUE)
    curses.init_pair(colors['cyan'], curses.COLOR_WHITE, curses.COLOR_CYAN)
    curses.init_pair(colors['grey'], curses.COLOR_WHITE, curses.COLOR_WHITE)
    curses.init_pair(colors['red'], curses.COLOR_WHITE, curses.COLOR_RED)

class Output:
    def __init__(self,name):
        self.name = name
        self.total = 0
        self.done = 0
        self.idle = 0
        self.t2 = 0
        self.t3 = 0
        self.missing = 0
    def add(self,state):
        self.total += 1
        setattr(self, state, getattr(self,state) + 1)
    def str(self):
        if self.total==0:
            return ''
        msgs = []
        s = ('%%-%is '%(WIDTH-1))%self.name[:(WIDTH-1)]
        if not args.monitor:
            msgs.append( s )
        else:
            msgs.append( (s,) )
        d_frac  = 1.*WIDTH*(self.done)/self.total
        t3_frac = 1.*WIDTH*(self.done+self.t3)/self.total
        t2_frac = 1.*WIDTH*(self.done+self.t3+self.t2)/self.total
        i_frac  = 1.*WIDTH*(self.done+self.t3+self.t2+self.idle)/self.total
        if not args.monitor:
            msgs.append( '[\033[0;%im'%colors['green'] )
        else:
            msgs.append( '[' )
        state = 0
        for i in xrange(WIDTH):
            if not args.monitor:
                s = ''
                if i>=i_frac:
                    if state != 4:
                        s += '\033[0;%im'%colors['red']
                        state = 4
                elif i>=t2_frac:
                    if state!=3:
                        s += '\033[0;%im'%colors['grey']
                        state = 3
                elif i>=t3_frac:
                    if state!=2:
                        s += '\033[0;%im'%colors['cyan']
                        state = 2
                elif i>=d_frac:
                    if state!=1:
                        s += '\033[0;%im'%colors['blue']
                        state = 1
                s += ' '
                msgs.append( s )
            else:
                if i>=i_frac:
                    msgs.append( (' ', curses.color_pair(colors['red'])) )
                elif i>=t2_frac:
                    msgs.append( (' ', curses.color_pair(colors['grey'])) )
                elif i>=t3_frac:
                    msgs.append( (' ', curses.color_pair(colors['cyan'])) )
                elif i>=d_frac:
                    msgs.append( (' ', curses.color_pair(colors['blue'])) )
                else:
                    msgs.append( (' ', curses.color_pair(colors['green'])) )
        s = ''
        if not args.monitor:
            s += '\033[0m] '
        else:
            s += ']'
        s += '%10i '%self.total
        s += '%10i '%self.t3
        s += '%10i '%self.t2
        s += '%10i '%self.idle
        s += '%10i '%self.missing
        s += '%10i '%self.done
        s += '(done=%.2f%%)\n'%(d_frac*100./WIDTH)
        if not args.monitor:
            msgs.append( s )
            return ''.join(msgs)
        else:
            msgs.append( (s,) )
            return msgs

def main(stdscr=None):
    if args.monitor:
        init_colors()
        stdscr.nodelay(True)
    while True:
        if args.monitor:
            c = stdscr.getch()
            curses.flushinp()
            if c == ord('q'):
                return

        global last_lock, last_check
        if time() - last_check > 5:
            if len(glob(lockdir+'/*')) > 0:
                recent_lock = int(path.getmtime(lockdir)) 
            else:
                recent_lock = 1
            last_check = time()

        if (recent_lock >= last_lock) or (time() - last_lock > args.monitor):
            # determine what files have been processed and logged as such
            processedfiles = []
            if args.verbose:
                print 'Finding locks...                      \r',
                sys.stdout.flush()
            locks = glob(lockdir+'/*lock')
            nl = len(locks)
            il = 1
            for lock in locks:
                if args.verbose:
                    print 'Reading lock %i/%i                   \r'%(il,nl),
                    sys.stdout.flush()
                il += 1
                flock = open(lock)
                for l in flock:
                    processedfiles.append(l.strip())

            if args.verbose:
                print 'Checking jobs...                 \r',
                sys.stdout.flush()

            # determine what samples from previous resubmissions are still running
            t2_samples = []
            t3_samples = []
            idle_samples = []
            if path.isfile(workdir+'/submission.pkl'): 
                with open(workdir+'/submission.pkl','rb') as fpkl:
                    submissions = pickle.load(fpkl)
            else:
                submissions = []
            for s in submissions:
                results = s.query_status()
                t3_samples += results['T3']
                t2_samples += results['T2']
                idle_samples += results['idle']

            t2_files = list(chain.from_iterable([x.files for x in t2_samples]))
            t3_files = list(chain.from_iterable([x.files for x in t3_samples]))
            idle_files = list(chain.from_iterable([x.files for x in idle_samples]))


            # for fancy display
            outputs = {}
            data = Output('Data')
            mc = Output('MC')

            if args.verbose:
                print 'Rebuilding configuration...            \r',

            all_samples = read_sample_config(args.infile)
            filtered_samples = {}
            merged_samples = {}
            outfile = open(args.outfile,'w')
            for name in sorted(all_samples):
                sample = all_samples[name]
                out_sample = DataSample(name,sample.dtype,sample.xsec)

                base_name = sub('_[0-9]+$','',name)
                if base_name not in outputs:
                    outputs[base_name] = Output(base_name)
                output = outputs[base_name]
                if base_name not in merged_samples:
                    merged_samples[base_name] = DataSample(base_name,sample.dtype,sample.xsec)
                merged_sample = merged_samples[base_name]

                to_resubmit = []

                for f in sample.files:
                    state = 'missing'
                    if f in processedfiles:
                        state = 'done'
                    elif f in t3_files:
                        state = 't3'
                    elif f in t2_files:
                        state = 't2'
                    elif f in idle_files:
                        state = 'idle'

                    if state=='missing' or (args.force and state!='done'):
                        out_sample.add_file(f)
                        merged_sample.add_file(f)

                    output.add(state)
                    if sample.dtype=='MC':
                        mc.add(state)
                    else:
                        data.add(state)

                if len(out_sample.files)>0:
                    filtered_samples[name] = out_sample

            if args.nfiles<0:
                keys = sorted(filtered_samples)
                for k in keys:
                    sample = filtered_samples[k]
                    if len(sample.files)==0:
                        continue
                    configs = sample.get_config(-1)
                    for c in configs:
                        outfile.write(c)
            else:
                keys = sorted(merged_samples)
                counter=0
                for k in keys:
                    sample = merged_samples[k]
                    if len(sample.files)==0:
                        continue
                    configs = sample.get_config(args.nfiles,suffix='_%i')
                    for c in configs:
                        outfile.write(c%(counter,counter))
                        counter += 1

            if args.verbose:
                print '\r',
                sys.stdout.flush()
            msg = []
            if args.verbose:
                msg.append('Summary:                     ')

            msg.append(header)
            if args.monitor:
                msg.append('\n')

            if args.monitor and len(outputs)+10>rows:
                args.silent = True
                msg.append( ('Too many samples to show in monitoring mode!\n', curses.color_pair(colors['red'])) )
            if not args.silent:
                for n in sorted(outputs):
                    if args.monitor:
                        msg.extend(outputs[n].str())
                    else:
                        msg.append(outputs[n].str().strip())
                msg.append('')
            if args.monitor:
                msg.extend(data.str())
                msg.extend(mc.str())
            else:
                msg.append(data.str().strip())
                msg.append(mc.str().strip())
            msg.append('')
            if args.monitor:
                msg.append('Legend: Done=[')
                msg.append( ('    ',curses.color_pair(colors['green'])) )
                msg.append('], T3=[')
                msg.append( ('    ',curses.color_pair(colors['blue'])) )
                msg.append('], T2=[')
                msg.append( ('    ',curses.color_pair(colors['cyan'])) )
                msg.append('], Idle[')
                msg.append( ('    ',curses.color_pair(colors['grey'])) )
                msg.append('], Missing=[')
                msg.append( ('    ',curses.color_pair(colors['red'])) )
                msg.append(']\n')
            else:
                msg.append( 'Legend: Done=\033[0;%im    \033[0m, T3=\033[0;%im    \033[0m, T2=\033[0;%im    \033[0m, Idle=\033[0;%im    \033[0m, Missing=\033[0;%im    \033[0m, '%(colors['green'],colors['blue'],colors['cyan'],colors['grey'],colors['red']))

            outfile.close()

            msg.append(strftime('%Y-%m-%d %H:%M:%S'))
            if args.monitor:
                msg.append('\n')
            msg.append( '\nMost recent submission:')
            if args.monitor:
                msg.extend([x+'\n' for x in query()])
                msg.append('\nPress q to close')
            else:
                msg.extend(query())
            msg.append('')

            if args.monitor:
                stdscr.clear()
                for m in msg:
                    if type(m) == str:
                        stdscr.addstr(m)
                    else:
                        stdscr.addstr(*m)
                stdscr.refresh()
            else:
                sys.stdout.write('\n'.join(msg))

            last_lock = int(time())

            if args.resubmit and len(merged_samples):
                resubmit()


        if args.monitor:
            sleep(.1)
        else:
            return


if args.monitor:
    curses.wrapper(main)
else:
    main()
