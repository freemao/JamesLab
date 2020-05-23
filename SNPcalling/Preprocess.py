    # -*- coding: UTF-8 -*-

"""
Prepare fastq files ready for SNP calling
"""
import re
import sys
import subprocess
import numpy as np
import pandas as pd
import os.path as op
from pathlib import Path
from subprocess import run
from schnablelab.apps.Tools import GenDataFrameFromPath
from schnablelab.apps.natsort import natsorted
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

def main():
    actions = (
        ('fastqc', 'check the reads quality'),
        ('trim_paired', 'quality control on paired reads'),
        ('trim_single', 'quality control on single reads'),
        ('combineFQ', 'combine splitted fastq files'),
        ('pre_ref', 'index the reference genome sequences'),
        ('pre_fqs', 'prepare fastq files read for mapping'),
        ('align_pe', 'paired-end alignment using bwa'),
        ('sam2bam', 'convert sam format to bam format'),
        ('sortbam', 'sort bam files'),
        ('index_bam', 'index bam files'),
        ('split_fa_region', 'genearte a list of freebayes/bamtools region specifiers'),
        ('bam_list', 'genearte a list of bam files for freebayes -L use'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def pre_ref(args):
    """
    %prog pre_ref ref.fa

    index the reference genome sequences using bwa, samtools, and picard tools
    """
    p = OptionParser(pre_ref.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=pre_ref.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref_fn, = args
    ref_fn, ref_dir = Path(ref_fn), Path(ref_fn).parent
    if not ref_fn.exists():
        sys.exit(f'reference file {ref_fn} does not exist!')
    ref_prefix = re.split('.fa|.fasta', ref_fn.name)[0]
    bwa_idx_exs = ('.amb', '.ann', '.bwt', '.pac', '.sa')
    bwa_bool = sum([(ref_dir/(ref_prefix+bie)).exists() for bie in bwa_idx_exs])
    cmds = []
    if bwa_bool !=5:
        print('bwa index does not exist...')
        cmd = f'ml bwa\nbwa index -p {ref_dir/ref_prefix} {ref_fn}'
        cmds.append(cmd)

    if not (ref_dir/(ref_fn.name+'.fai')).exists():
        print('fai index does not exist...')
        cmd = f'ml samtools\nsamtools faidx {ref_fn}'
        cmds.append(cmd)
    
    dict_fn = ref_dir/(ref_prefix+'.dict')
    if not dict_fn.exists():
        print('dict index does not exist...')
        cmd = f'ml gatk4\ngatk CreateSequenceDictionary -R {ref_fn} -O {dict_fn}'
        cmds.append(cmd)

    if len(cmds)>0:
        if not opts.disable_slurm:
            put2slurm_dict = vars(opts)
            put2slurm(cmds, put2slurm_dict)
        else:
            print('commands running on local:\n%s'%('\n'.join(cmds)))
    else:
        print('All reference index files have already existed!')

def find_sm(target_str, re_pattern):
    sms = re_pattern.findall(target_str)
    if len(sms)==1:
        sm = sms[0][1:-1]
        return '-'.join(re.split('[_-]', sm))
    else:
        sys.exit(f"bad file name '{target_str}'!")

def pre_fqs(args):
    """
    %prog pre_fqs dir1 dir2 ... output.csv

    parse RG and SM info of all fastq files and get them ready for mapping

    dir1: where fastq files are located
        add more directories if fastq files are located at different directories
    output.csv:
        output csv file containing all parsed fq files
    """
    p = OptionParser(pre_fqs.__doc__)
    p.add_option('--fq_fn_pattern', default='*.fastq.gz',
                help = 'file extension of fastq files')
    p.add_option('--sm_re_pattern', default=r"[^a-z0-9]P[0-9]{3}[_-]W[A-Z][0-9]{2}[^a-z0-9]", 
                help = 'the regular expression pattern to pull sample name from filename')
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    *fq_dirs, out_csv = args

    tmp_df_ls = []
    for fq_dir in fq_dirs:
        fq_dir = Path(fq_dir)
        if not fq_dir.exists():
            sys.exit(f'{fq_dir} does not exist!')
        tmp_df = GenDataFrameFromPath(fq_dir, pattern=opts.fq_fn_pattern)
        if tmp_df.shape[0]==0:
            sys.exit(f"no fastq files found under '{fq_dir}' directory!")
        print(f'{tmp_df.shape[0]} fastq files found under {fq_dir}!')
        tmp_df_ls.append(tmp_df)
    df = pd.concat(tmp_df_ls)

    prog = re.compile(opts.sm_re_pattern)
    df['sm'] = df['fn'].apply(lambda x: find_sm(x, prog))
    df = df.sort_values(['sm', 'fn']).reset_index(drop=True)
    print(f"Total {df['sm'].unique().shape[0]} samples found!")
    print(df['sm'].value_counts())
    df.to_csv(out_csv, index=False)
    print(f'{out_csv} has been generated')

def align_pe(args):
    """
    %prog align_pe ref_indx_base fq_fns.csv

    paire-end alignment using bwa.
    args:
        ref_index_base: the prefix of reference index files
        fq_fns.csv: the csv file including parsed fq files from pre_fqs function.
    """
    p = OptionParser(align_pe.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=align_pe.__name__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    ref_base, fq_csv = args
    df = pd.read_csv(fq_csv)

    df_R1, df_R2 = df[::2], df[1::2]
    if df_R1.shape[0] != df_R2.shape[0]:
        sys.exit('number of R1 and R2 files are not consistent!')

    cmds = []
    for (_,r1), (_,r2) in zip(df_R1.iterrows(), df_R2.iterrows()):
        r1_fn, r2_fn, sm = r1['fnpath'], r2['fnpath'], r1['sm']
        r1_fn_arr, r2_fn_arr = np.array(list(r1_fn.name)), np.array(list(r2_fn.name))
        bools = (r1_fn_arr != r2_fn_arr)
        if bools.sum() != 1:
            print(r1_fn, r2_fn)
            sys.exit('check fq file names!')
        idx = np.argmax(bools)
        prefix = re.split('[-_]R', r1_fn.name[:idx])[0]
        RG = r"'@RG\tID:%s\tSM:%s'"%(sm, sm)
        bam_fn = f'{prefix}.pe.sorted.bam'
        cmd = f"bwa mem -t {opts.ncpus_per_node} -R {RG} {ref_base} {r1_fn} {r2_fn} | samtools sort -@{opts.ncpus_per_node} -o {bam_fn} -"
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml bwa\nml samtools'
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['cmd_header'] = cmd_header
        put2slurm(cmds, put2slurm_dict)

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
        header = Slurm_header%(10, 8000, prf, prf, prf)
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
