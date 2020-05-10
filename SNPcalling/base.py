# -*- coding: UTF-8 -*-

"""
base class and functions to handle with vcf file
"""
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from schnablelab.apps.base import ActionDispatcher, OptionParser

def main():
    actions = (
        ('FilterMissing', 'filter out SNPs with high missing rate'),
        ('FilterMAF', 'filter out SNP with extremely low minor allele frequency'),
        ('FilterHetero', 'filter out SNPs with high heterozygous rates'),
        ('SubsamplingSNPs', 'grep a subset of SNPs from a vcf file'),
        ('SubsamplingSMs', 'grep a subset of samples from a vcf file'),
        ('SummarizeLD', 'ld decay in log scale'),
        ('vcf2hmp', 'convert vcf to hmp foramt'),
        ('Info', 'get basic information of a vcf file')
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
            num_SNPs = 0
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
                    num_SNPs += 1
        self.num_SNPs = num_SNPs
        self.SMs_header = SMs_header
        self.SMs = SMs
        self.numSMs = numSMs
        self.numHash = n
        self.HashChunk = hash_chunk
        self.HashChunk2 = hash_chunk2
        self.numHeaderLines = len(self.HashChunk)
        self.hmpfield11 = 'rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode'
        self.hmpheader = self.hmpfield11 + '\t' + '\t'.join(self.SMs) + '\n'

    def AsDataframe(self):
        df = pd.read_csv(self.fn, skiprows=range(self.numHash-1), delim_whitespace=True)
        return df
    
    def ToHmp(self):
        '''
        yield line in hmp format
        '''
        cen_NA = '+\tNA\tNA\tNA\tNA\tNA\tNA'
        with open(self.fn) as f:
            for _ in range(self.numHash):
                next(f)
            for i in f:
                j = i.split()
                
                a1, a2 = j[3], j[4]
                if len(a1) == len(a2) ==1:
                    a1a2 = ''.join([a1, a2])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = a1*2, a2*2
                elif len(a1)==1 and len(a2)>1:
                    a1a2 = ''.join(['-', a2[-1]])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = '--', a2[-1]*2
                elif len(a1) >1 and len(a2)==1:
                    a1a2 = ''.join([a1[-1], '-'])
                    a2a1 = a1a2[::-1]
                    a1a1, a2a2 = a1[-1]*2, '--'
                else:
                    print('bad format line:\n  %s'%j[2])
                    continue
                geno_dict = {'0/0':a1a1, '0|0':a1a1, 
                            '0/1':a1a2, '0|1':a1a2, 
                            '1/0':a2a1, '1|0':a2a1,
                            '1/1':a2a2, '1|1':a2a2,
                            './.':'NN', '.|.':'NN'}
                genos = list(map(geno_dict.get, j[9:]))
                if None in genos:
                    print(i)
                    sys.exit('unknow genotype detected!')
                genos = '\t'.join(genos)
                rs, chr, pos = j[2], j[0], j[1]
                alleles = '/'.join([a1a2[0], a1a2[1]])
                new_line = '\t'.join([rs, alleles, chr, pos, cen_NA, genos])+'\n'
                yield new_line

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

def vcf2hmp(args):
    '''
    %prog vcf2hmp input_vcf
    convert file in vcf format to hmp format
    can also try: 'run_pipeline.pl -Xms512m -Xmx10G -fork1 -vcf vcf_fn -export -exportType Hapmap'
    '''
    p = OptionParser(vcf2hmp.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputhmp = Path(inputvcf).name.replace('.vcf', '.hmp')

    vcf = ParseVCF(inputvcf)
    with open(outputhmp, 'w') as f:
        f.write(vcf.hmpheader)
        pbar = tqdm(vcf.ToHmp(), total=vcf.num_SNPs, desc='vcf 2 hmp', position=0)
        for i in pbar:
            f.write(i)
            pbar.set_description('converting %s'%i.split()[2])
            #pbar.update(1)
    print('Done! check output %s...'%outputhmp)    

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
    outputvcf = Path(inputvcf).name.replace('.vcf', '_mis%s.vcf'%opts.missing_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.Missing(), total=vcf.num_SNPs, desc='Filter Missing', position=0)
        for i, miss in pbar:
            if miss <= opts.missing_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing %s'%i.split()[0])
    print('Done! %s SNPs removed! check output %s...'%(n, outputvcf))

def FilterMAF(args):
    """
    %prog FilterMAF input_vcf
    Remove rare MAF SNPs
    """
    p = OptionParser(FilterMAF.__doc__)
    p.add_option('--maf_cutoff', default = 0.01, type='float',
        help = 'specify the MAF rate cutoff, SNPs lower than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    outputvcf = Path(inputvcf).name.replace('.vcf', '_maf%s.vcf'%opts.maf_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.MAF(), total=vcf.num_SNPs, desc='Filter MAF', position=0)
        for i, maf in pbar:
            if maf >= opts.maf_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing %s'%i.split()[0])
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
    outputvcf = Path(inputvcf).name.replace('.vcf', '_het%s.vcf'%opts.het_cutoff)

    vcf = ParseVCF(inputvcf)
    n = 0
    with open(outputvcf, 'w') as f:
        f.writelines(vcf.HashChunk)
        pbar = tqdm(vcf.Hetero(), total=vcf.num_SNPs, desc='Filter Heterozygous', position=0)
        for i, het in pbar:
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing %s'%i.split()[0])
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
    outputvcf = Path(inputvcf).name.replace('.vcf', '_subSNPs.vcf')

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
    outputvcf = Path(inputvcf).name.replace('.vcf', '_subSMs.vcf')

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

def Info(args):
    """
    %prog Info input_vcf

    get basic info including SMs, number of SMs, number of SNPs, 
    """
    p = OptionParser(Info.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputvcf, = args
    vcf = ParseVCF(inputvcf)
    print('number of samples: {val:,}'.format(val=vcf.numSMs))
    print("number of hash ('#') lines: {val:,}".format(val=vcf.numHeaderLines))
    print('number of SNPs: {val:,}'.format(val=vcf.num_SNPs))
    print('Sample names: \n  %s'%'\n  '.join(vcf.SMs))

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
