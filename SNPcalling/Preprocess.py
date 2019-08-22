# -*- coding: UTF-8 -*-

"""
Prepare fastq files ready for SNP calling
"""

import sys
import os.path as op
from pathlib import Path
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob
from schnablelab.apps.natsort import natsorted
import subprocess
from schnablelab.apps.headers import Slurm_header


def main():
    actions = (
        ('fastqc', 'check the reads quality'),
        ('trim_paired', 'quality control on paired reads'),
        ('trim_single', 'quality control on single reads'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def fastqc(args):
    """
    %prog fastqc in_dir out_dir
        in_dir: the dir where fastq files are located
        out_dir: the dir saving fastqc reports

    generate slurm files for fastqc jobs
    """
    p = OptionParser(fastqc.__doc__)
    p.add_option("--pattern", default = '*.fastq', 
            help="the pattern of fastq files, qutation needed") 
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args

    out_path = Path(out_dir)
    if not out_path.exists():
        print('%s does not exist...')
    dir_path = Path(in_dir)
    fqs = dir_path.glob(opts.pattern)
    for fq in fqs:
        prf = '.'.join(fq.name.split('.')[0:-1])
        print(prf)
        cmd = 'fastqc %s -o %s'%(str(fq), out_dir)
        header = Slurm_header%(10, 10000, prf, prf, prf)
        header += 'ml fastqc\n'
        header += cmd
        with open('%s.fastqc.slurm'%(prf), 'w') as f:
            f.write(header)

def trim_paired(args):
    """
    %prog trim in_dir out_dir
    quality control on the paired reads
    """
    p = OptionParser(trim_paired.__doc__)
    p.add_option('--pattern_r1', default = '*_R1.fastq',
            help='filename pattern for forward reads')
    p.add_option('--pattern_r2', default = '*_R2.fastq',
            help='filename pattern for reverse reads')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir,out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        print('output dir %s does not exist...'%out_dir)
        sys.exit(1)
    r1_fns = glob('%s/%s'%(in_dir, opts.pattern_r1))
    r2_fns = glob('%s/%s'%(in_dir, opts.pattern_r2))
    for r1_fn, r2_fn in zip(r1_fns, r2_fns):
        r1_path = Path(r1_fn)
        r2_path = Path(r2_fn)
        prf = '_'.join(r1_path.name.split('_')[0:-1])+'.PE'
        print(prf)
        r1_fn_out1 = r1_path.name.replace('R1.fastq', 'trimed.R1.fastq')
        r1_fn_out2 = r1_path.name.replace('R1.fastq', 'unpaired.R1.fastq')
        r2_fn_out1 = r2_path.name.replace('R2.fastq', 'trimed.R2.fastq')
        r2_fn_out2 = r2_path.name.replace('R2.fastq', 'unpaired.R2.fastq')
        cmd = 'java -jar $TM_HOME/trimmomatic.jar PE -phred33 %s %s %s %s %s %s TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40'%(r1_fn,r2_fn,str(out_path/r1_fn_out1),str(out_path/r1_fn_out2),str(out_path/r2_fn_out1),str(out_path/r2_fn_out2))
        header = Slurm_header%(10, 10000, prf, prf, prf)
        header += 'ml trimmomatic\n'
        header += cmd
        with open('%s.trim.slurm'%(prf), 'w') as f:
            f.write(header)

def trim_single(args):
    """
    %prog trim in_dir out_dir
    quality control on the single end reads
    """
    p = OptionParser(trim_paired.__doc__)
    p.add_option('--pattern', default = '*_Unpaired.fastq',
            help='filename pattern for all single end reads')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir,out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        print('output dir %s does not exist...'%out_dir)
        sys.exit(1)
    fns = glob('%s/%s'%(in_dir, opts.pattern))
    for fn in fns:
        fn_path = Path(fn)
        prf = '_'.join(fn_path.name.split('_')[0:-1])+'.SE'
        print(prf)
        fn_out = fn_path.name.replace('Unpaired.fastq', 'trimed.Unpaired.fastq')
        cmd = 'java -jar $TM_HOME/trimmomatic.jar SE -phred33 %s %s TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40'%(fn, str(out_path/fn_out))
        header = Slurm_header%(10, 10000, prf, prf, prf)
        header += 'ml trimmomatic\n'
        header += cmd
        with open('%s.trim.slurm'%(prf), 'w') as f:
            f.write(header)

if __name__ == "__main__":
    main()
