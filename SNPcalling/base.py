# -*- coding: UTF-8 -*-

"""
Filter SNPs using bcftools.
Find more details at bcftools website:
https://samtools.github.io/bcftools/bcftools.html
"""

import pandas as pd
import numpy as np
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
        ('Missing', 'filter missing rate using customized script'),
        ('MAF', 'filter minor allele frequency using customized script'),
        ('Heterozygous', 'filter SNPs with high heterozygous rates'),
        ('Bad_Indels', 'remove wrong INDELs'),
        ('GrepImputatedVcf', 'grep the SNPs with lower missing rate before imputation from whole imputed vcf'),
        ('SummarizeLD', 'ld decay in log scale'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def SummarizeLD(args):
    """
    %prog ld.csv num0 out.txt
    ld.csv: ld tab delimited file generated from tassel
    num0: 0s in the distance

    summarize ld decay in log scale 0-100kb
    """
    p = OptionParser(SummarizeLD.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ld_fn,num0,out_fn, = args
    df = pd.read_csv(ld_fn, delim_whitespace=True, usecols=['Dist_bp', 'R^2'])
    df = df.dropna().sort_values('Dist_bp').reset_index(drop=True)

    mybin = [10**i for i in np.arange(0, float(num0)+0.1, 0.1)]
    blockPreIndex = np.histogram(df['Dist_bp'].values, bins=mybin)[0]

    a = list(blockPreIndex)
    a.insert(0,0)
    boxlist = []
    for idx,ele in enumerate(a):
        st = sum(a[0:idx])
        ed = sum(a[0:idx+1])
        boxlist.append(df['R^2'][st:ed].values)
    boxlist.pop(0)
    
    with open(out_fn, 'w') as f:
        for idx,ele in enumerate(boxlist):
            if len(ele) >= 1:
                averageR2, sd = sum(ele)/float(len(ele)), np.var(ele)
            elif len(ele) == 0:
                averageR2, sd = '',''
            f.write('%s\t%s\t%s\t%s\n'%(10**(idx*0.1),(10**((idx+1)*0.1)), averageR2, sd))

def Missing(args):
    """
    %prog vcf
    Remove SNPs with high missing rate
    """
    p = OptionParser(Missing.__doc__)
    p.add_option('--missing_rate', default = 0.7, type='float', 
        help = 'specify the missing rate cutoff. SNPs with missing rate higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, = args
    prefix = vcffile.split('.vcf')[0]
    new_f = prefix + '.mis%s.vcf'%opts.missing_rate
    total_n = getSMsNum(vcffile)
    print('Total %s Samples.'%total_n)
    with open(new_f, 'w') as f_out:
        with open(vcffile) as f_in:
            for i in f_in:
                if i.startswith('#'):
                    f_out.write(i)
                else:
                    _,_,_,m = genotypes_count(i)
                    missing_rate = m/total_n
                    if missing_rate <= opts.missing_rate:
                        f_out.write(i)

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
    %prog vcf_in
    Remove bad and high heterizygous loci
    """
    p = OptionParser(Heterozygous.__doc__)
    p.add_option('--h2_rate', default = 0.1, type='float',
        help = 'specify the heterozygous rate cutoff, higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, = args
    prefix = vcffile.split('.vcf')[0]
    new_f = prefix + '.hete%s.vcf'%opts.h2_rate
    f0 = open(vcffile)
    f1 = open(new_f, 'w')
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
    call('ml bcftools', shell=True)
    child = subprocess.Popen('bcftools query -l %s|wc -l'%vcffile, shell=True, stdout=subprocess.PIPE)
    SMs_num = int(child.communicate()[0])
    return SMs_num

def MAF(args):
    """
    %prog MAF vcf maf

    filter rare MAF SNPs
    """
    p = OptionParser(MAF.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, maf, = args
    prefix = vcffile.split('.vcf')[0]
    new_f = prefix + '.maf%s.vcf'%maf

    total_n = getSMsNum(vcffile)
    print('Total %s Samples.'%total_n)

    with open(new_f, 'w') as f_out:
        with open(vcffile) as f_in:
            for i in f_in:
                if i.startswith('#'):
                    f_out.write(i)
                else:
                    ref, alt, het, mis = genotypes_count(i)
                    an1, an2 = ref*2+het, alt*2+het
                    maf = min(an1,an2)/(an1+an2)
                    if maf >= float(maf):
                        f_out.write(i)
    
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
