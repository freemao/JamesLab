# -*- coding: UTF-8 -*-

"""
Prepare fastq files ready for SNP calling
"""

import sys
import subprocess
import os.path as op
from pathlib import Path
from subprocess import run
from schnablelab.apps.natsort import natsorted
from schnablelab.apps.headers import Slurm_header
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob


def main():
    actions = (
        ('fastqc', 'check the reads quality'),
        ('trim_paired', 'quality control on paired reads'),
        ('trim_single', 'quality control on single reads'),
        ('combineFQ', 'combine splitted fastq files'),
        ('index_ref', 'index the genome sequences'),
        ('align', 'align reads to the reference'),
        ('sam2bam', 'convert sam format to bam format'),
        ('sortbam', 'sort bam files'),
        ('index_bam', 'index bam files'),
        ('split_fa_region', 'genearte a list of freebayes/bamtools region specifiers'),
        ('bam_list', 'genearte a list of bam files for freebayes -L use'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bam_list(args):
    """
    %prog bam_list bam_dir out_fn

    genearte a list of bam files for freebayes -L use
    """
    p = OptionParser(bam_list.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    bam_dir, fn_out, = args
    dir_path = Path(bam_dir)
    bams = sorted(dir_path.glob('*.bam'))
    f = open(fn_out, 'w')
    for bam in bams:
        f.write('%s\n'%bam)
    f.close()

def split_fa_region(args):
    """
    %prog fa.fai region_size out_fn
        fa.fai: index file for the fa file
        region_size: the size for each splitted region
        out_fn: the output file

    genearte a list of freebayes/bamtools region specifiers
    """
    p = OptionParser(split_fa_region.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    fasta_index_file, region_size, fn_out, = args
    fasta_index_file = open(fasta_index_file)
    region_size = int(region_size)
    fn_out = open(fn_out, 'w')
    for line in fasta_index_file:
        fields = line.strip().split("\t")
        chrom_name = fields[0]
        chrom_length = int(fields[1])
        region_start = 0
        while region_start < chrom_length:
            start = region_start
            end = region_start + region_size
            if end > chrom_length:
                end = chrom_length
            line = chrom_name + ":" + str(region_start) + "-" + str(end)+'\n'
            fn_out.write(line)
            region_start = end
    fn_out.close()

def index_bam(args):
    """
    %prog bam_dir 
        bam_dir: sorted bam files folder

    index bam files using samtools/0.1 sort function.
    """
    p = OptionParser(index_bam.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    bam_dir, = args
    dir_path = Path(bam_dir)
    bams = dir_path.glob('*.sorted.bam')
    for bam in bams:
        prf = bam.name.split('.sorted.bam')[0]
        cmd = 'samtools index %s'%bam
        header = Slurm_header%(10, 10000, prf, prf, prf)
        header += 'ml samtools/0.1\n'
        header += cmd
        with open('%s.indexbam.slurm'%prf, 'w') as f:
            f.write(header)

def sortbam(args):
    """
    %prog in_dir out_dir
        in_dir: bam files folder
        out_dir: sorted bam files folder

    sort bam files using samtools/0.1 sort function.
    """
    p = OptionParser(sortbam.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args

    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    bams = dir_path.glob('*.bam')
    for bam in bams:
        prf = bam.name.split('.bam')[0]
        sort_bam = prf+'.sorted'
        sort_bam_path = out_path/sort_bam
        cmd = 'samtools sort %s %s'%(bam, sort_bam_path)
        header = Slurm_header%(100, 15000, prf, prf, prf)
        header += 'ml samtools/0.1\n'
        header += cmd
        with open('%s.sortbam.slurm'%prf, 'w') as f:
            f.write(header)

def sam2bam(args):
    """
    %prog in_dir out_dir
        in_dir: sam files folder
        out_dir: bam files folder

    convert sam to bam using samtools/0.1.
    """
    p = OptionParser(sam2bam.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args

    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    sams = dir_path.glob('*.sam')
    for sam in sams:
        prf = sam.name.split('.sam')[0]
        bam = prf+'.bam'
        bam_path = out_path/bam
        cmd = 'samtools view -bS %s > %s'%(sam, bam_path)
        header = Slurm_header%(100, 15000, prf, prf, prf)
        header += 'ml samtools/0.1\n'
        header += cmd
        with open('%s.sam2bam.slurm'%prf, 'w') as f:
            f.write(header)

def align(args):
    """
    %prog align indx_base fq_fn ...

    do alignment using bwa.
    """
    p = OptionParser(align.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    ref_base = args[0]
    fq_fns = args[1:]
    print(fq_fns)
    sm = Path(fq_fns[0]).name.split('_trim')[0]
    gid = sm.split('R')[0]
    print(gid)
    R = r"'@RG\tID:%s\tSM:%s'"%(gid, gid)
    if len(fq_fns)==1:
        sam = sm+'.se.sam'
        print('run single-end alignment')
        cmd = 'bwa mem -R %s %s %s > %s \n'%(R, ref_base, fq_fns[0], sam)
        prf = '%s.se.align'%sm
    elif len(fq_fns)==2:
        sam = sm+'.pe.sam'
        print('run paired-end alignment')
        cmd = 'bwa mem -R %s %s %s %s > %s \n'%(R, ref_base, fq_fns[0], fq_fns[1], sam)
        prf = '%s.pe.align'%sm
    else:
        sys.exit('only one or two read files')
    header = Slurm_header%(100, 10000, prf, prf, prf)
    header += 'ml bwa\n'
    header += cmd
    with open('%s.slurm'%prf, 'w') as f:
        f.write(header)

def index_ref(args):
    """
    %prog index_ref ref.fa

    index the reference genome sequences
    """
    p = OptionParser(index_ref.__doc__)
    p.add_option('--tool', default='bwa', choices=('bwa', 'samtools'),
            help = 'tool for indexing reference genome')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref_fn, = args
    prefix = '.'.join(ref_fn.split('.')[0:-1])
    if opts.tool == 'bwa':
        cmd = 'bwa index -p %s %s'%(prefix, ref_fn)
        print(cmd)
        header = Slurm_header%(100, 15000, prefix, prefix, prefix)
        header += 'ml bwa\n'
        header += cmd
        with open('%s.bwa_index.slurm'%prefix, 'w') as f:
            f.write(header)
    else:
        cmd = 'samtools faidx %s'%ref_fn
        print(cmd)
        header = Slurm_header%(10, 10000, prefix, prefix, prefix)
        header += 'ml samtools\n'
        header += cmd
        with open('%s.samtools_index.slurm'%prefix, 'w') as f:
            f.write(header)

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
        sys.exit('%s does not exist...')
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
        sys.exit('output dir %s does not exist...'%out_dir)
    r1_fns = glob('%s/%s'%(in_dir, opts.pattern_r1))
    r2_fns = glob('%s/%s'%(in_dir, opts.pattern_r2))
    for r1_fn, r2_fn in zip(r1_fns, r2_fns):
        r1_path = Path(r1_fn)
        r2_path = Path(r2_fn)
        prf = '_'.join(r1_path.name.split('_')[0:-1])+'.PE'
        print(prf)
        r1_fn_out1 = r1_path.name.replace('R1.fastq', 'trim.R1.fastq')
        r1_fn_out2 = r1_path.name.replace('R1.fastq', 'unpaired.R1.fastq')
        r2_fn_out1 = r2_path.name.replace('R2.fastq', 'trim.R2.fastq')
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
        sys.exit('output dir %s does not exist...'%out_dir)
    fns = glob('%s/%s'%(in_dir, opts.pattern))
    for fn in fns:
        fn_path = Path(fn)
        prf = '_'.join(fn_path.name.split('_')[0:-1])+'.SE'
        print(prf)
        fn_out = fn_path.name.replace('Unpaired.fastq', 'trim.Unpaired.fastq')
        cmd = 'java -jar $TM_HOME/trimmomatic.jar SE -phred33 %s %s TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40'%(fn, str(out_path/fn_out))
        header = Slurm_header%(10, 10000, prf, prf, prf)
        header += 'ml trimmomatic\n'
        header += cmd
        with open('%s.trim.slurm'%(prf), 'w') as f:
            f.write(header)

def combineFQ(args):
    """
    %prog combineFQ pattern(with quotation) fn_out
    """

    p = OptionParser(combineFQ.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    fq_pattern, fn_out, = args
    fns = glob(fq_pattern)
    cmd = 'cat %s > %s'%(' '.join(fns), fn_out)
    print(cmd)
    run(cmd, shell=True)

if __name__ == "__main__":
    main()
