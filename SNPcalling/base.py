# -*- coding: UTF-8 -*-

"""
base class and functions to handle with vcf file
"""
import sys
import numpy as np
import pandas as pd
from schnablelab.apps.base import ActionDispatcher, OptionParser

def main():
    actions = (
        ('FilterMissing', 'filter out SNPs with high missing rate'),
        ('FilterMAF', 'filter out SNP with extremely low minor allele frequency'),
        ('FilterHetero', 'filter out SNPs with high heterozygous rates'),
        ('SubsamplingSNPs', 'grep a subset of SNPs from a vcf file'),
        ('SubsamplingSMs', 'grep a subset of samples from a vcf file'),
        ('SummarizeLD', 'ld decay in log scale'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())

class ParseVCF():
    '''
    parse vcf file
    '''
    def __init__(self, filename):
        self.fn = filename
        with open(filename) as f:
            n = 0
            hash_chunk = []
            hash_chunk2 = []
            for i in f:
                if i.startswith('##'):
                    n += 1
                    hash_chunk.append(i)
                    hash_chunk2.append(i)
                    continue
                if i.startswith('#'):
                    SMs_header = i.split()[:9]
                    SMs = i.split()[9:]
                    numSMs = len(SMs)
                    n += 1
                    hash_chunk.append(i)
                else:
                    break
        self.SMs_header = SMs_header
        self.SMs = SMs
        self.numSMs = numSMs
        self.numHash = n
        self.HashChunk = hash_chunk
        self.HashChunk2 = hash_chunk2

    def AsDataframe(self):
        df = pd.read_csv(self.fn, skiprows=range(self.numHash-1), delim_whitespace=True)
        return df

    def Missing(self):
        '''
        yield missing rate for each line
        '''
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_miss = i.count('./.')+i.count('.|.')
                yield i, num_miss/self.numSMs

    def Hetero(self):
        '''
        yield (line, heterozygous rate) for each line
        '''
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_a = i.count('0/0')+i.count('0|0')
                num_b = i.count('1/1')+ i.count('1|1')
                num_h = i.count('0/1')+ i.count('0|1')+i.count('1|0')
                if num_h > max(num_a, num_b):
                    yield i, 1
                else:
                    yield i, num_h/float(num_a + num_b + num_h)

    def MAF(self):
        '''
        yield minor allele frequence for each line
        '''
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                num_a = i.count('0/0')+i.count('0|0')
                num_b = i.count('1/1')+ i.count('1|1')
                num_h = i.count('0/1')+ i.count('0|1')+i.count('1|0')
                a1, a2 = num_a*2+num_h, num_b*2+num_h
                yield  i, min(a1,a2)/(a1+a2)

def FilterMissing(args):
    """
    %prog FilterMissing input_vcf
    Remove SNPs with high missing rate
    """
    p = OptionParser(FilterMissing.__doc__)
    p.add_option('--missing_cutoff', default = 0.7, type='float', 
        help = 'specify the missing rate cutoff. SNPs higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = inputvcf.split('.vcf')[0] + '_mis%s.vcf'%opts.missing_cutoff

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        for i, miss in vcf.Missing():
            if miss <= opts.missing_cutoff:
                f.write(i)
            else:
                n += 1
    print('Done! %s SNPs removed! check output %s...'%(n, outputvcf))

def FilterHetero(args):
    """
    %prog FilterHetero input_vcf
    Remove bad and high heterizygous loci
    """
    p = OptionParser(FilterHetero.__doc__)
    p.add_option('--het_cutoff', default = 0.1, type='float',
        help = 'specify the heterozygous rate cutoff, SNPs higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = inputvcf.split('.vcf')[0] + '_het%s.vcf'%opts.het_cutoff

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        for i, het in vcf.Hetero():
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
    print('Done! %s SNPs removed! check output %s...'%(n, outputvcf))

def FilterMAF(args):
    """
    %prog FilterMAF input_vcf
    Remove rare MAF SNPs
    """
    p = OptionParser(FilterMAF.__doc__)
    p.add_option('--MAF_cutoff', default = 0.01, type='float',
        help = 'specify the MAF rate cutoff, SNPs lower than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = inputvcf.split('.vcf')[0] + '_maf%s.vcf'%opts.MAF_cutoff

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        for i, maf in vcf.MAF():
            if maf >= opts.MAF_cutoff:
                f.write(i)
            else:
                n += 1
    print('Done! %s SNPs removed! check output %s...'%(n, outputvcf))
    
def SubsamplingSNPs(args):
    """
    %prog SubsamplingSNPs input_vcf SNPs.csv 
    grep a subset of SNPs defined in SNPs.csv (One ID per row without header) from the input_vcf
    """
    p = OptionParser(SubsamplingSNPs.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, SNPcsv, = args
    outputvcf = inputvcf.split('.vcf')[0] + '_subSNPs.vcf'

    vcf = ParseVCF(inputvcf)
    df_vcf = vcf.AsDataframe()

    IDs = pd.read_csv(SNPcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print('number of specified SNPs: %s'%num_IDs)
    df_vcf = df_vcf[df_vcf['ID'].isin(IDs)]
    print('%s out of %s found in VCF'%(df_vcf.shape[0], num_IDs))
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk2)
    df_vcf.to_csv(outputvcf, sep='\t', index=False, mode='a')
    print('Done! check output %s...'%outputvcf)

def SubsamplingSMs(args):
    """
    %prog SubsamplingSMs input_vcf SMs.csv 
    grep a subset of samples defined in SMs.csv (One sample name per row without header) from the input_vcf
    """
    p = OptionParser(SubsamplingSMs.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, SMcsv, = args
    outputvcf = inputvcf.split('.vcf')[0] + '_subSMs.vcf'

    vcf = ParseVCF(inputvcf)
    df_vcf = vcf.AsDataframe()

    IDs = pd.read_csv(SMcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print('number of specified Samples: %s'%num_IDs)

    subsm = vcf.SMs_header
    for id in IDs:
        if id not in vcf.SMs:
            print('%s not found in vcf...'%id)
        else:
            subsm.append(id)
    print('%s out of %s found in VCF'%(len(subsm)-9, num_IDs))

    df_vcf = df_vcf[subsm]
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk2)
    df_vcf.to_csv(outputvcf, sep='\t', index=False, mode='a')
    print('Done! check output %s...'%outputvcf)   

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

if __name__ == "__main__":
    main()
