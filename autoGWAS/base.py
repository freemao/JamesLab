# -*- coding: UTF-8 -*-

"""
base class and functions to handle with hmp file and GWAS results
"""
import re
import sys
import numpy as np
import pandas as pd
import os.path as op
from tqdm import tqdm
from pathlib import Path
from subprocess import call
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

plink = op.abspath(op.dirname(__file__)) + '/../apps/plink'
GEC = op.abspath(op.dirname(__file__)) + '/../apps/gec.jar'

def main():
    actions = (
        ('FilterMissing', 'filter out SNPs with high missing rate'),
        ('FilterMAF', 'filter out SNP with extremely low minor allele frequency'),
        ('FilterHetero', 'filter out SNPs with high heterozygous rates'),
        ('SubsamplingSNPs', 'grep a subset of specified SNPs from a hmp file'),
        ('DownsamplingSNPs', 'pick up some SNPs from a huge hmp file using Linux sed command'),
        ('SubsamplingSMs', 'grep a subset of samples from a hmp file'),
        ('hmp2ped', 'convert hmp file to plink map and ped file'),
        ('ped2bed', 'convert plink ped format to binary bed format'),
        ('IndePvalue', 'estimate the number of independent SNPs using GEC'),
        ('HmpSingle2Double', 'convert single hmp to double type hmp'),
        ('Info', 'get basic info for a hmp file'),
        ('MAFs', 'calculate the MAF for all SNPs in hmp'),
        ('SortHmp', 'Sort hmp file based on chromosome order and position')
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())

# N:missing, -:gap
geno_codification = {'A':'AA', 'C':'CC', 'G':'GG', 'T':'TT',
    'R':'AG', 'Y':'CT', 'S':'GC', 'W':'AT', 'K':'GT', 'M':'AC',
    'N':'NN', '-':'--'} 

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
            #print('guess hmp type: %s'%type)
            numSNPs = sum(1 for _ in f)
        self.headerline = headerline
        self.SMs_header = SMs_header
        self.SMs = SMs
        self.numSMs = numSMs
        self.numSNPs = numSNPs+1
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
        if self.type=='single':
            print('converting the single type to double type...')
            df_part1 = df.iloc[:, :11]
            df_part2 = df.iloc[:, 11:].applymap(geno_codification.get).fillna('NN')
            df = pd.concat([df_part1, df_part2], axis=1)
            self.type = 'double'
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
        pbar = tqdm(self.SMs)
        for col in pbar:
            col_a1, col_a2 = df_hmp[col].str.get(0), df_hmp[col].str.get(1)
            col_a1.index = np.arange(0, df_hmp.shape[0]*2, 2)
            col_a2.index = np.arange(1, df_hmp.shape[0]*2, 2)
            tmp_s= pd.concat([col_a1, col_a2]).sort_index()
            tmp_ss.append(tmp_s)
            pbar.set_description('converting %s'%col)
        df_ped_part2 = pd.DataFrame(tmp_ss).reset_index(drop=True)

        df_ped = pd.concat([df_ped_part1, df_ped_part2], axis=1)
        return df_map, df_ped
    
    def BIMBAM(self):
        pass
    
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
                        try:
                            maf = min(a1, a2)/(a1+a2)
                        except ZeroDivisionError:
                            yield i, 0
                        else:
                            yield i, maf 
        else:
            with open(self.fn) as f:
                next(f)
                NN_ls = ['NN' for i in range(self.numSMs)]
                for i in f:
                    j = i.split()
                    alleles = j[1].split('/')
                    try:
                        allele1, allele2 = alleles
                    except ValueError:
                        yield i, 0
                    else:
                        genos_single = j[11:]
                        genos_double = list(map(geno_codification.get, genos_single, NN_ls))
                        genos = ''.join(genos_double)
                        a1, a2 = genos.count(allele1), genos.count(allele2)
                        try:
                            maf = min(a1, a2)/(a1+a2)
                        except ZeroDivisionError:
                            yield i, 0
                        else:
                            yield i, maf 

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
    outputhmp = Path(inputhmp).name.replace('.hmp', '_mis%s.hmp'%opts.missing_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.Missing(), total=hmp.numSNPs)
        for i, miss in pbar:
            if miss <= opts.missing_cutoff:
                f.write(i)
            else: 
                n +=1
            pbar.set_description('processing %s'%i.split()[2])
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
    outputhmp = Path(inputhmp).name.replace('.hmp', '_het%s.hmp'%opts.het_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.Hetero(), total=hmp.numSNPs)
        for i, het in pbar:
            if het <= opts.het_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing %s'%i.split()[2])
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
    outputhmp = Path(inputhmp).name.replace('.hmp', '_maf%s.hmp'%opts.MAF_cutoff)

    hmp = ParseHmp(inputhmp)
    n = 0
    with open(outputhmp, 'w') as f:
        f.write(hmp.headerline)
        pbar = tqdm(hmp.MAF(), total=hmp.numSNPs)
        for i, maf in pbar:
            if maf >= opts.MAF_cutoff:
                f.write(i)
            else:
                n += 1
            pbar.set_description('processing %s'%i.split()[2])
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
    outputhmp = Path(inputhmp).name.replace('.hmp', '_subSNPs.hmp')

    hmp = ParseHmp(inputhmp)
    df_hmp = hmp.AsDataframe()

    IDs = pd.read_csv(SNPcsv, header=None)[0].values
    num_IDs = IDs.shape[0]
    print('number of specified SNPs: %s'%num_IDs)
    df_hmp = df_hmp[df_hmp['rs#'].isin(IDs)]
    print('%s out of %s found in Hmp'%(df_hmp.shape[0], num_IDs))
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
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
    outputhmp = Path(inputhmp).name.replace('.hmp', '_subSMs.hmp')

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
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print('Done! check output %s...'%outputhmp)

def DownsamplingSNPs(args):
    """
    %prog downsampling input_hmp

    Pick up some SNPs from a huge hmp file using Linux sed command
    """
    p = OptionParser(DownsamplingSNPs.__doc__)
    p.add_option('--downscale', default=10,
                 help='specify the downscale level')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=DownsamplingSNPs.__name__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_ds%s.hmp'% opts.downsize)
    cmd = "sed -n '1~%sp' %s > %s" % (opts.downsize, inputhmp, outputhmp)
    print('cmd:\n%s\n' % cmd)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm([cmd], put2slurm_dict)

def hmp2ped(args):
    """
    %prog input_hmp

    Convert hmp file to Plink map and ped files
    """
    p = OptionParser(hmp2ped.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    output_prefix = Path(inputhmp).name.rstrip('.hmp')

    hmp = ParseHmp(inputhmp)
    df_map, df_ped = hmp.AsMapPed(missing=False)
    df_map.to_csv('%s.map'%output_prefix, sep='\t', index=False, header=None)
    df_ped.to_csv('%s.ped'%output_prefix, sep='\t', index=False, header=None)

def ped2bed(args):
    """
    %prog ped_prefix

    Convert plink ped/map to binary bed/bim/fam format using Plink
    """
    p = OptionParser(ped2bed.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='add this option to disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=ped2bed.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ped_prefix, = args
    cmd_header = 'ml plink'
    cmd = 'plink --noweb --file %s --make-bed --out %s' % (ped_prefix, ped_prefix)
    print('cmd on HCC:\n%s\n%s' % (cmd_header, cmd))

    cmd_local = '%s --noweb --file %s --make-bed --out %s' % (plink, ped_prefix, ped_prefix)
    print('cmd on local desktop:\n%s\n'%cmd_local)
    
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['cmd_header'] = cmd_header
        put2slurm([cmd], put2slurm_dict)

def IndePvalue(args):
    """
    %prog IndePvalue bed_prefix output_fn

    Estimate number of idenpendent SNPs using GEC
    """
    p = OptionParser(IndePvalue.__doc__)
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='add this option to disable converting commands to slurm jobs')
    p.add_slurm_opts(job_prefix=IndePvalue.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    bed_prefix, output_fn = args
    cmd = 'java -Xmx18g -jar %s --noweb --effect-number --plink-binary %s --genome --out %s' % (GEC, bed_prefix, output_fn)
    print('cmd:\n%s\n' % cmd)

    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['memory'] = 20000
        put2slurm([cmd], put2slurm_dict)

def HmpSingle2Double(args):
    """
    %prog HmpSingle2Double input_single_hmp 
    convert single type hmp file to double type hmp file
    """
    p = OptionParser(HmpSingle2Double.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_db.hmp')

    hmp = ParseHmp(inputhmp)
    df_hmp = hmp.AsDataframe()
    df_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print('Done! check output %s...'%outputhmp)

def Info(args):
    """
    %prog Info input_hmp
    get basic info for a hmp file
    """
    p = OptionParser(Info.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    hmp = ParseHmp(inputhmp)

    print('Genotype type: %s'%hmp.type)
    print('Number of samples: {val:,}'.format(val=hmp.numSMs))
    print('Number of SNPs: {val:,}'.format(val=hmp.numSNPs))
    print('Sample names: \n  %s'%'\n  '.join(hmp.SMs))

def MAFs(args):
    """
    %prog MAFs input_hmp 

    calculate MAF for all SNPs in hmp
    """
    p = OptionParser(MAFs.__doc__)
    _, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputcsv = Path(inputhmp).name.replace('.hmp', '.maf.csv')
    hmp = ParseHmp(inputhmp)
    with open(outputcsv, 'w') as f:
        pbar = tqdm(hmp.MAF(), total=hmp.numSNPs, desc='get MAF', position=0)
        for i, maf in pbar:
            f.write('%s\n'%maf)
            pbar.set_description('calculating %s'%i.split()[2])
    print('Done! check output %s...'%(outputcsv))

def SortHmp(args):
    """
    %prog SortHmp input_hmp 
    Sort hmp based on chromosome order and position. Can also try tassel:
     'run_pipeline.pl -Xms16g -Xmx18g -SortGenotypeFilePlugin -inputFile in_fn -outputFile out_fn -fileType Hapmap'
    """
    p = OptionParser(SortHmp.__doc__)
    _, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    inputhmp, = args
    outputhmp = Path(inputhmp).name.replace('.hmp', '_sorted.hmp')

    hmp = ParseHmp(inputhmp)
    df_sorted_hmp = hmp.AsDataframe(needsort=True)
    df_sorted_hmp.to_csv(outputhmp, sep='\t', index=False, na_rep='NA')
    print('Done! check output %s...'%outputhmp)

if __name__ == "__main__":
    main()