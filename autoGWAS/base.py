# -*- coding: UTF-8 -*-

"""
base class and functions used in GWAS results analysis
"""
import re
import sys
import numpy as np
import pandas as pd

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
    def __init__(self, filename, type='double'):
        '''
        args:
            filename: hmp file name
            type: hmp format. double or single
        '''
        self.fn = filename
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
    
    def MAFs(self):
        '''
        calculate MAF for each SNP
        return: MAFs in pandas series format
        '''
        allmafs = []
        if self.type == 'double':
            with open(self.fn) as f:
                f.readline()
                for i in f:
                    j = i.split()
                    print(j[0])
                    alleles = j[1].split('/')
                    if len(alleles) == 2:
                        allele1, allele2 = alleles
                        genos = ''.join(j[11:])
                        a1, a2 = genos.count(allele1), genos.count(allele2)
                        maf = a1/float(a1+a2) \
                            if a1 <= a2 \
                            else a2/float(a1+a2)
                        allmafs.append(maf)
                    else:
                        allmafs.append(0)
        else:
            sys.exit('not implement yet..')
        allmafs = pd.Series(allmafs)
        return allmafs

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

