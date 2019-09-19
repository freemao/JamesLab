# -*- coding: UTF-8 -*-

"""
Filter SNPs using bcftools.
Find more details at bcftools website:
https://samtools.github.io/bcftools/bcftools.html
"""

import pandas as pd
import os.path as op
import sys
from pathlib import Path
from subprocess import call
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob,iglob
from schnablelab.apps.natsort import natsorted
import subprocess
from schnablelab.apps.headers import Slurm_header

def main():
    actions = (
        ('NUM_ALT', 'filter number of alternative SNPs'),
        ('Missing', 'filter missing rate'),
        ('Heterozygous', 'filter SNPs with high heterozygous rates'),
        ('Bad_Indels', 'remove wrong INDELs'),
        ('MAF_unimputed', 'filter minor allele frequency using bcftools'),
        ('Subsampling', 'choose part of samples from vcf'),
        ('GrepImputatedVcf', 'grep the SNPs with lower missing rate before imputation from whole imputed vcf'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def Subsampling(args):
    """
    %prog Subsampling SMs_file vcf_or_vcf.gz
    Subsampling vcf file using bcftools. The samples order will also change following the order in SMs_file.
    """
    p = OptionParser(Subsampling.__doc__)
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    SMsfile, vcffile, = args
    
    prefix = vcffile.split('/')[-1].split('.vcf')[0]
    new_f = prefix + '.subsm.vcf'
    cmd = "bcftools view -S %s %s > %s\n"%(SMsfile, vcffile, new_f)
    print(cmd)
    jobfile = '%s.subsm.slurm'%prefix
    f = open(jobfile, 'w')
    header = Slurm_header%(opts.time, opts.memory, opts.prefix, opts.prefix, opts.prefix)
    header += 'module load bcftools\n'
    header += cmd
    f.write(header) 
    print('slurm file %s.subsm.slurm has been created, you can sbatch your job file.'%prefix)

def NUM_ALT(args):
    """
    %prog NUM_ALT vcf_or_vcf.gz
    only retain SNPs with only one ALT    
    """
    p = OptionParser(NUM_ALT.__doc__)
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)
    
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, = args
    prefix = vcffile.split('.')[0]
    new_f = prefix + '.alt1.vcf'
    cmd = "bcftools view -i 'N_ALT=1' %s > %s"%(vcffile, new_f)
    jobfile = '%s.alt1.slurm'%prefix
    f = open(jobfile, 'w')
    header = Slurm_header%(opts.time, opts.memory, opts.prefix, opts.prefix, opts.prefix)
    header += 'module load bacftools\n'
    header += cmd
    f.write(header) 
    print('slurm file %s.alt1.slurm has been created, you can sbatch your job file.'%prefix)

def Missing(args):
    """
    %prog vcf_or_vcf.gz
    Remove SNPs with high missing rate (>0.4 by default)
    """
    p = OptionParser(Missing.__doc__)
    p.add_option('--missing_rate', default = 0.4, 
        help = 'specify the missing rate cutoff. SNPs with missing rate higher than this cutoff will be removed.')
    p.add_option('--NS', default = 'NS', 
        help = 'specify the tag name to calculate the number of nonmissing samples. If NS, NZ are unavilable, can specify AN cause AN/2==NS')
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, = args
    prefix = vcffile.split('.vcf')[0]
    new_f = prefix + '.mis.vcf'
    nonMR = 1-float(opts.missing_rate)

    out = getSMsNum(vcffile)
    print('Total %s Samples.'%out)

    if opts.NS in ('NS', 'NZ'):
        cmd = "bcftools view -i '%s/%s >= %.2f' %s > %s\n"%(opts.NS, out, nonMR, vcffile, new_f)
    elif opts.NS == 'AN':
        cmd = "bcftools view -i 'AN/%s >= %.2f' %s > %s\n"%(out*2, nonMR, vcffile, new_f)
    else:
        sys.exit('NS, NZ, AN only')

    jobfile = '%s.mis.slurm'%prefix
    f = open(jobfile, 'w')
    header = Slurm_header%(opts.time, opts.memory, opts.prefix, opts.prefix, opts.prefix)
    header += 'module load bcftools\n'
    header += cmd
    f.write(header)
    print('slurm file %s.missR.slurm has been created, you can sbatch your job file.'%prefix)

def genotypes_count(snp):
    """
    calculate the numbers of ref, alt, hetero, missing genotypes.
    """
    a1 = snp.count('0/0')+snp.count('0|0')
    a2 = snp.count('1/1')+ snp.count('1|1')
    h = snp.count('0/1')+ snp.count('0|1')+snp.count('1|0')
    m = snp.count('./.')+snp.count('.|.')
    return a1, a2, h, m

def Heterozygous(args):
    """
    %prog vcf_in vcf_out
    Remove bad and high heterizygous loci
    """
    p = OptionParser(Heterozygous.__doc__)
    p.add_option('--h2_rate', default = 0.05,
        help = 'specify the heterozygous rate cutoff, higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, vcf_out, = args
    prefix = vcffile.split('.vcf')[0]
    f0 = open(vcffile)
    f1 = open(vcf_out, 'w')
    for i in f0:
        if i.startswith('#'):
            f1.write(i)
        else:
            a1, a2, h, m = genotypes_count(i)
            if h <= max(a1, a2) and h/float(a1+a2+h) <= float(opts.h2_rate):
                f1.write(i)
    f0.close()
    f1.close()

def getSMsNum(vcffile):
    call('module load bcftools', shell=True)
    child = subprocess.Popen('bcftools query -l %s|wc -l'%vcffile, shell=True, stdout=subprocess.PIPE)
    SMs_num = float(child.communicate()[0])
    return SMs_num

def MAF_unimputed(args):
    """
    %prog MAF_unimputed in_dir out_dir
    do filtering on MAF using bcftools.
    """
    p = OptionParser(MAF_unimputed.__doc__)
    p.add_option("--pattern", default = '*.vcf', 
            help="the pattern of vcf files, qutation needed") 
    p.add_option('--maf', default = '0.01',
        help = 'specify the missing rate cutoff, MAF smaller than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args   
    maf = float(opts.maf)

    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        prf = vcf.name.split('.vcf')[0]
        new_vcf = prf +'.maf.vcf'
        cmd = "bcftools view -i 'AF>=%s && AF<=%s' %s > %s"%(maf, 1-maf, vcf, out_path/new_vcf) 
        header = Slurm_header%(10, 5000, prf, prf, prf)
        header += 'ml bcftools\n'
        header += cmd
        with open('%s.maf%s.slurm'%(prf, maf), 'w') as f:
            f.write(header)

def MAF_imputed(args):
    '''
    else:
        f1 = open(new_f, 'w')
        with open(vcffile) as f0:
            for i in f0:
                if i.startswith('#'): f1.write(i)
                else:
                    ref, alt, het, mis = genotypes_count(i)
                    an1, an2 = ref*2+het, alt*2+het
                    maf = min(an1,an2)/(an1+an2)
                    if maf >= float(opts.maf):
                        f1.write(i)
        f1.close()
    '''
    
def GrepImputatedVcf(args):
    """
    %prog LowerMissingVcf ImputedVcf out_vcf
    grep SNPs with lower missing rate before imputation from whole imputed vcf
    """
    p = OptionParser(GrepImputatedVcf.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    low_vcf, ipt_vcf,out_vcf = args

    seed_head_n = 0
    with open(low_vcf) as f:
        for i in f:
            if i.startswith('##'): seed_head_n += 1
            else: break
    low_vcf_df = pd.read_csv(low_vcf, delim_whitespace=True, usecols=['#CHROM', 'POS'], skiprows=seed_head_n)
    seed = low_vcf_df['#CHROM']+'\t'+low_vcf_df['POS'].astype('str')
    print('seed generated.')

    ipt_head_n = 0
    with open(ipt_vcf) as f1:
        for i in f1:
            if i.startswith('##'): ipt_head_n += 1
            else: break
    ipt_vcf_df = pd.read_csv(ipt_vcf, delim_whitespace=True, skiprows=ipt_head_n)
    target = ipt_vcf_df['#CHROM']+'\t'+ipt_vcf_df['POS'].astype('str')
    print('whole imputed target generated.')

    seed_bool = target.isin(seed)
    out_vcf_df = ipt_vcf_df[seed_bool]
    out_vcf_df.to_csv(out_vcf, index=False, sep='\t')
    print('Done! check %s'%out_vcf)

if __name__ == "__main__":
    main()
