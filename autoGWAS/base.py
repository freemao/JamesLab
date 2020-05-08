# -*- coding: UTF-8 -*-

"""
base class and functions to handle with hmp file and GWAS results
"""
import re
import sys
import numpy as np
import pandas as pd
from subprocess import call
from schnablelab.apps.base import ActionDispatcher, OptionParser

def main():
    actions = (
        ('FilterMissing', 'filter out SNPs with high missing rate'),
        ('FilterMAF', 'filter out SNP with extremely low minor allele frequency'),
        ('FilterHetero', 'filter out SNPs with high heterozygous rates'),
        ('SubsamplingSNPs', 'grep a subset of specified SNPs from a hmp file'),
        ('downsamplingSNPs', 'grep a subset of SNPs from a hmp file'),
        ('SubsamplingSMs', 'grep a subset of samples from a hmp file'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())

geno_codification = {
    'A':'AA',
    'C':'CC',
    'G':'GG',
    'T':'TT',
    'R':'AG',
    'Y':'CT',
    'S':'GC',
    'W':'AT',
    'K':'GT',
    'M':'AC',
    'N':'NN'}

def sortchr(x):
    '''
    criteria to sort chromosome names
    '''
    x1 = re.findall(r'\d+$', x)
    if len(x1)==1:
        return int(x1[0])
    else:
        sys.exit('check chromosome name!')

class ParseHmp():
    '''
    parse hmp file
    '''
    def __init__(self, filename):
        '''
        args:
            filename: hmp file name
            type: hmp format. double or single
        '''
        self.fn = filename
        with open(filename) as f:
            headerline = f.readline()
            SMs_header = headerline.split()[:11]
            SMs = headerline.split()[11:]
            numSMs = len(SMs)
            firstgeno = f.readline().split()[11]
            type = 'single' if len(firstgeno)==1 else 'double'
            print('guess hmp type: %s'%type)
        self.headerline = headerline
        self.SMs_header = SMs_header
        self.SMs = SMs
        self.numSMs = numSMs
        self.type = type
        
    def AsDataframe(self, needsort=False):
        '''
        args:
            needsort: if hmp need to be sorted
        '''
        df = pd.read_csv(self.fn, delim_whitespace=True)
        df['chrom'] = df['chrom'].astype('str')
        if needsort:
            chrs = list(df['chrom'].unique())
            chrs_ordered = sorted(chrs, key=sortchr)
            df['chrom'] = pd.Categorical(df['chrom'], chrs_ordered, ordered=True)
            df = df.sort_values(['chrom', 'pos']).reset_index(drop=True)
        return df
    
    def AsMapPed(self, missing=False):
        df_hmp = self.AsDataframe()
        df_map = df_hmp[['rs#', 'chrom', 'pos']]
        df_map['centi'] = 0
        map_cols = ['chrom', 'rs#', 'centi', 'pos']
        df_map = df_map[map_cols]

        part1_cols = ['fam_id', 'indi_id', 'indi_id_father', 'indi_id_mother', 'sex', 'pheno']
        zeros = np.zeros(self.numSMs, dtype=int)
        df_ped_part1 = pd.DataFrame(dict(zip(part1_cols, [zeros, self.SMs, zeros, zeros, zeros, zeros])) )
        
        df_hmp = df_hmp.iloc[:, 11:]
        if missing:
            df_hmp = df_hmp.replace('NN', '00')
        tmp_ss = []
        for col in df_hmp:
            col_a1, col_a2 = df_hmp[col].str.get(0), df_hmp[col].str.get(1)
            col1.index = np.arange(0, df_hmp.shape[0]*2, 2)
            col2.index = np.arange(1, df_hmp.shape[0]*2, 2)
            tmp_s= pd.concat([col_a1, col_a2]).sort_index()
            tmp_ss.append(tmp_s)
        df_ped_part2 = pd.DataFrame(tmp_ss).reset_index(drop=True)

        df_ped = pd.concat([df_ped_part1, df_ped_part2])
        return df_map, df_ped
    
    def Missing(self):
        '''
        yield (line, missing rate) for each line
        '''
        if self.type == 'double':
            with open(self.fn) as f:
                next(f)
                for i in f:
                    num_miss = i.split()[11:].count('NN')
                    yield i, num_miss/self.numSMs
        else:
            with open(self.fn) as f:
                next(f)
                for i in f:
                    num_miss = i.split()[11:].count('N')
                    yield i, num_miss/self.numSMs

    def MAF(self):
        '''
        yield (line, maf) for each line
        '''
        if self.type == 'double':
            with open(self.fn) as f:
                next(f)
                for i in f:
                    j = i.split()
                    alleles = j[1].split('/')
                    try:
                        allele1, allele2 = alleles
                    except ValueError:
                        yield i, 0
                    else:
                        genos = ''.join(j[11:])
                        a1, a2 = genos.count(allele1), genos.count(allele2)
                        yield i, min(a1, a2)/(a1+a2)
        else:
            with open(self.fn) as f:
                next(f)
                for i in f:
                    j = i.split()
                    alleles = j[1].split('/')
                    try:
                        allele1, allele2 = alleles
                    except ValueError:
                        yield i, 0
                    else:
                        genos = j[11:]
                        genos = list(map(geno_codification.get, genos, genos))
                        if None in genos:
                            print(i)
                            sys.exit('unkown character in hmp file')
                        genos = ''.join(genos)
                        a1, a2 = genos.count(allele1), genos.count(allele2)
                        yield i, min(a1, a2)/(a1+a2)

    def Hetero(self):
        '''
        yield (line, heterozgous rate) for each line
        '''
        if self.type == 'double':
            with open(self.fn) as f:
                next(f)
                for i in f:
                    j = i.split()
                    genos = j[11:]
                    alleles = j[1].split('/')
                    ab, ba = ''.join(alleles), ''.join(list(reversed(alleles)))
                    aa, bb = list(map(lambda x: x*2, alleles))
                    num_a = genos.count(aa)
                    num_b = genos.count(bb)
                    num_h = genos.count(ab)+genos.count(ba)
                    if num_h > max(num_a, num_b):
                        yield i, 1
                    else:
                        yield i, num_h/float(num_a + num_b + num_h)
        else:
            with open(self.fn) as f:
                next(f)
                for i in f:
                    j = i.split()
                    genos = j[11:]
                    genos = list(map(geno_codification.get, genos, genos))
                    if None in genos:
                            print(i)
                            sys.exit('unkown character in hmp file')
                    alleles = j[1].split('/')
                    ab, ba = ''.join(alleles), ''.join(list(reversed(alleles)))
                    aa, bb = list(map(lambda x: x*2, alleles))
                    num_a = genos.count(aa)
                    num_b = genos.count(bb)
                    num_h = genos.count(ab)+genos.count(ba)
                    if num_h > max(num_a, num_b):
                        yield i, 1
                    else:
                        yield i, num_h/float(num_a + num_b + num_h)

class ReadGWASfile():
    
    def __init__(self, filename, software, needsort=False, usecols=None):
        '''
        Args:
            filename: gwas result filename
            software: gwas software (gemma, farmcpu, mvp, gapit)
            needsort: if the gwas file need to be sorted
            mvp_p_col: specify which pvlue column if multiple approaches used in MVP
        '''
        self.fn = filename
        self.software = software
        self.needsort = needsort
        self.usecols = usecols
        
        if self.software == 'gemma':
            df = pd.read_csv(self.fn, delim_whitespace=True, usecols=['chr', 'rs', 'ps', 'p_lrt'])
            df = df[['rs', 'chr', 'ps', 'p_lrt']]
        elif self.software == 'farmcpu':
            df = pd.read_csv(self.fn, usecols=['SNP', 'Chromosome', 'Position', 'P.value'])
        elif self.software == 'gapit':
            df = pd.read_csv(self.fn, usecols=['SNP', 'Chromosome', 'Position ', 'P.value'])
        elif self.software == 'other':
            if self.usecols is None:
                sys.exit('specify which columns for use if choosing other!')
            if (not isinstance(self.usecols, list)):
                sys.exit('usecols must be a list')
            if len(self.usecols) != 4:
                sys.exit('usecols must have the lenght of 4')
            df = pd.read_csv(self.fn, usecols=self.usecols)
        else:
            sys.exit('only gemma, farmcpu, gapit, and other are supported!')
        df.columns = ['snp', 'chr', 'pos', 'pvalue']
        df['chr'] = df['chr'].astype('str')
        df['pvalue'] = -np.log10(df['pvalue'])
        df.columns = ['snp', 'chr', 'pos', '-log10Pvalue']
        if needsort:
            chrs = list(df['chr'].unique())
            chrs_ordered = sorted(chrs, key=sortchr)
            df['chr'] = pd.Categorical(df['chr'], chrs_ordered, ordered=True)
            df = df.sort_values(['chr', 'pos']).reset_index(drop=True)
        
        self.df = df
        self.numberofSNPs = df.shape[0]

    def SignificantSNPs(self, p_cutoff=0.05, MeRatio=1):
        '''
        extract Significant SNPs
        '''
        cutoff = -np.log10(p_cutoff/(MeRatio * self.numberofSNPs))
        df_sigs = self.df[self.df['-log10Pvalue'] >= cutoff].reset_index(drop=True)
        return df_sigs

def FilterMissing(args):
    """
    %prog FilterMissing input_hmp
    Remove SNPs with high missing rate
    """
    p = OptionParser(FilterMissing.__doc__)
    p.add_option('--missing_cutoff', default = 0.7, type='float', 
        help = 'specify the missing rate cutoff. SNPs higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = inputhmp.replace('.hmp', '_mis%s.hmp'%opts.missing_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        for i, miss in hmp.Missing():
            if miss <= opts.missing_cutoff:
                f.write(i)
            else: 
                n +=1
    print('Done! %s SNPs removed! check output %s...'%(n, outputhmp))

def FilterHetero(args):
    """
    %prog FilterHetero input_hmp
    Remove bad and high heterizygous loci (coducting Missing and MAF first)
    """
    p = OptionParser(FilterHetero.__doc__)
    p.add_option('--het_cutoff', default = 0.1, type='float',
        help = 'specify the heterozygous rate cutoff, SNPs higher than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = inputhmp.replace('.hmp', '_het%s.hmp'%opts.het_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.writelines(hmp.headerline)
        for i, het in hmp.Hetero():
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
    print('Done! %s SNPs removed! check output %s...'%(n, outputhmp))

def FilterMAF(args):
    """
    %prog FilterMAF input_hmp
    Remove rare MAF SNPs (conducting Missing filter first)
    """
    p = OptionParser(FilterMAF.__doc__)
    p.add_option('--MAF_cutoff', default = 0.01, type='float',
        help = 'specify the MAF rate cutoff, SNPs lower than this cutoff will be removed.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = inputhmp.replace('.hmp', '_maf%s.hmp'%opts.MAF_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.writelines(hmp.headerline)
        for i, maf in hmp.MAF():
            if maf >= opts.MAF_cutoff:
                f.write(i)
            else:
                n += 1
    print('Done! %s SNPs removed! check output %s...'%(n, outputhmp))
    
def SubsamplingSNPs(args):
    """
    %prog SubsamplingSNPs input_hmp SNPs.csv 
    grep a subset of SNPs defined in SNPs.csv (One ID per row without header) from the input_hmp
    """
    p = OptionParser(SubsamplingSNPs.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, SNPcsv, = args
    outputhmp = inputhmp.replace('.hmp', '_subSNPs.hmp')

    hmp = ParseHmp(inputhmp)
    df_hmp = hmp.AsDataframe()

    IDs = pd.read_csv(SNPcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print('number of specified SNPs: %s'%num_IDs)
    df_hmp = df_hmp[df_hmp['rs#'].isin(IDs)]
    print('%s out of %s found in Hmp'%(df_hmp.shape[0], num_IDs))
    df_hmp.to_csv(outputhmp, sep='\t', index=False)
    print('Done! check output %s...'%outputhmp)

def SubsamplingSMs(args):
    """
    %prog SubsamplingSMs input_hmp SMs.csv 
    grep a subset of samples defined in SMs.csv (One sample name per row without header) from the input_hmp
    """
    p = OptionParser(SubsamplingSMs.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, SMcsv, = args
    outputhmp = inputhmp.replace('.hmp', '_subSMs.hmp')

    hmp = ParseHmp(inputhmp)
    df_hmp = hmp.AsDataframe()

    IDs = pd.read_csv(SMcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print('number of specified Samples: %s'%num_IDs)

    subsm = hmp.SMs_header
    for id in IDs:
        if id not in hmp.SMs:
            print('%s not found in hmp...'%id)
        else:
            subsm.append(id)
    print('%s out of %s found in Hmp'%(len(subsm)-11, num_IDs))

    df_hmp = df_hmp[subsm]
    df_hmp.to_csv(outputhmp, sep='\t', index=False)
    print('Done! check output %s...'%outputhmp)

def downsamplingSNPs(args):
    """
    %prog downsampling input_hmp

    Choose part of SNPs as mapping markers when the genotype dataset is huge
    """
    p = OptionParser(downsamplingSNPs.__doc__)
    p.add_option('--downscale', default=10,
                 help='specify the downscale level')
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    inputhmp, = args
    outputhmp = inputhmp.replace('.hmp', '.ds%s.hmp'% opts.downsize)
    cmd = "sed -n '1~%sp' %s > %s" % (opts.downsize, inputhmp, outputhmp)
    call(cmd, shell=True)

if __name__ == "__main__":
    main()