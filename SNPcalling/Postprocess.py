# -*- coding: UTF-8 -*-

"""
Split a big file using sed.
Find more details at cnblog:
www.cnblogs.com/freemao/p/7076127.html
"""

from pathlib import Path
import os.path as op
import sys
import pandas as pd
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob
from schnablelab.apps.natsort import natsorted
import subprocess
from subprocess import run
from schnablelab.apps.headers import Slurm_header, multiCPU_header

# the location of linkimpute, beagle executable
lkipt = op.abspath(op.dirname(__file__)) + '/../apps/LinkImpute.jar'
begle = op.abspath(op.dirname(__file__)) + '/../apps/beagle.24Aug19.3e8.jar'
tassel = op.abspath(op.dirname(__file__)) + '/../apps/tassel-5-standalone/run_pipeline.pl'


def main():
    actions = (
        ('splitVCF', 'split a vcf to several smaller files with equal size'),
        ('merge_files', 'combine split vcf or hmp files'),
        ('combineFQ', 'combine split fqs'),
        ('impute_beagle', 'impute vcf using beagle or linkimpute'),
        ('vcf2hmp', 'convert vcf to hmp format'),
        ('FixIndelHmp', 'fix the indels problems in hmp file converted from tassel'),
        ('FilterVCF', 'remove bad snps using bcftools'),
        ('onlyMissing', 'perform filtering only on missing data'),
        ('rmHetero', 'remove high heterozygous loci'),
        ('subsampling', 'subsampling and reorder vcf files'),
        ('fixGTsep', 'fix the allele separator for beagle imputation'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def fixGTsep(args):
    """
    %prog fixGTsep in_dir out_dir

    replace the allele separator . in freebayes vcf file to / which is required for beagle
    """
    p = OptionParser(fixGTsep.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.fixGT.vcf'
        out_fn_path = out_path/out_fn
        cmd = "perl -pe 's/\s\.:/\t.\/.:/g' %s > %s"%(vcf, out_fn_path)
        header = Slurm_header%(10, 10000, sm, sm, sm)
        header += cmd
        with open('%s.fixGT.slurm'%sm, 'w') as f:
            f.write(header)

def subsampling(args):
    """
    %prog subsampling in_dir out_dir samples.csv

    subsampling and reorder samples in vcf file using bcftools
    """
    p = OptionParser(subsampling.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir,sm_fn, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.subsm.vcf'
        out_fn_path = out_path/out_fn
        cmd = 'bcftools view -S %s %s > %s'%(sm_fn, vcf, out_fn_path)
        header = Slurm_header%(10, 20000, sm, sm, sm)
        header += 'ml bcftools\n'
        header += cmd
        with open('%s.subsm.slurm'%sm, 'w') as f:
            f.write(header)

def rmHetero(args):
    """
    %prog rmHetero in_dir out_dir

    filter SNPs with high heterozygous rate
    """

    p = OptionParser(rmHetero.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_option('--rate', default='0.05',
                 help='heterozygous rate cutoff')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.rmHete.vcf'
        cmd = 'python -m schnablelab.SNPcalling.FilterSNPs Heterozygous %s %s --h2_rate %s'%(vcf, out_path/out_fn, opts.rate)
        header = Slurm_header%(10, 20000, sm, sm, sm)
        #header += 'conda activate MCY\n'
        header += cmd
        with open('%s.rmHete.slurm'%sm, 'w') as f:
            f.write(header)

def FilterVCF(args):
    """
    %prog FilterVCF dir_in dir_out

    filter SNPs using bcftools
    """

    p = OptionParser(FilterVCF.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_option('--n_alt', default='1',
                 help='number of alt')
    p.add_option('--qual',
                 help='minimum snp quality, 10: 10% is wrong, 20: 1% is wrong')
    p.add_option('--maf',
                 help='cutoff of minor allele frequency')
    p.add_option('--missing',
                 help='cutoff of missing rate')
    p.add_option('--stype', default='snps,indels',
                 help='snp types, comma separated if multiple types specified')
    p.add_option('--normalization', default=False, action='store_true',
                 help='perform normalization')
    p.add_option('--ref', default='/work/schnablelab/cmiao/TimeSeriesGWAS/Genotype_GBS/Reference_Genome_4th/Sbicolor_454_v3.0.1.fa',
                 help='required if --normalization specified')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)

    cond1 = 'N_ALT==%s'%opts.n_alt
    if opts.qual:
        cond1 += ' && QUAL>=%s'%opts.qual
    if opts.maf:
        cond1 += ' && MAF>=%s'%opts.maf
    if opts.missing:
        missing_rate = 1 - float(opts.missing)
        cond1 += ' && NS/N_SAMPLES > %.2f'%missing_rate
    cmd = "bcftools view -i '{cond1}' -v '{stype}' %s".format(cond1=cond1, stype=opts.stype)
    if opts.normalization:
        cmd += ' | bcftools norm -f %s -m -both'%(opts.ref)
    cmd += ' > %s'
    
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.bcflt.vcf'
        out_fn_path = out_path/out_fn
        header = Slurm_header%(10, 8000, sm, sm, sm)
        header += 'ml bcftools\n'
        header += cmd%(vcf, out_fn_path)
        with open('%s.bcflt.slurm'%sm, 'w') as f:
            f.write(header)

def getSMsNum(vcffile):
    call('module load bcftools', shell=True)
    child = subprocess.Popen('bcftools query -l %s|wc -l'%vcffile, shell=True, stdout=subprocess.PIPE)
    SMs_num = float(child.communicate()[0])
    return SMs_num


def onlyMissing(args):
    """
    %prog dir_in dir_out
    Remove SNPs with high missing rate (>0.4 by default)
    """
    p = OptionParser(onlyMissing.__doc__)
    p.add_option('--pattern', default = '*.vcf',
        help = 'remove loci with high missing data rate')
    p.add_option('--missing_rate', default = 0.7,
        help = 'specify the missing rate cutoff. SNPs with missing rate higher than this cutoff will be removed.')
    p.add_option('--software', default = 'GATK',
        help = 'specify the software generating the vcf file.')
    p.add_option('--total_sm', default = 'infer',
        help = 'specify total samples in the vcf file.')

    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args

    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    nonMR = 1-float(opts.missing_rate)

    for vcf in vcfs:
        prefix = '.'.join(vcf.name.split('.')[0:-1])
        total_sm = getSMsNum(str(vcf)) if opts.total_sm == 'infer' else opts.total_sm
        print('Total %s Samples.'%total_sm)
        out_fn = prefix+'.mis.vcf'
        out_fn_path = out_path/out_fn
        if opts.software == 'GATK':
            cmd = "bcftools view -i 'AN/%s >= %.2f' %s > %s\n"%(int(total_sm)*2, nonMR, vcf, out_fn_path)
        else:
            cmd = "bcftools view -i 'NS/%s >= %.2f' %s > %s\n"%(int(total_sm)*2, nonMR, vcf, out_fn_path)
        header = Slurm_header%(10, 8000, prefix, prefix, prefix)
        header += 'ml bcftools\n'
        header += cmd
        with open('%s.mis.slurm'%prefix, 'w') as f:
            f.write(header)

def splitVCF(args):
    """
    %prog splitVCF N vcf
    split vcf to N smaller files with equal size
    """
    p = OptionParser(splitVCF.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    N, vcffile, = args
    N = int(N)
    prefix = vcffile.split('.')[0]
    cmd_header = "sed -ne '/^#/p' %s > %s.header" % (vcffile, prefix)
    subprocess.call(cmd_header, shell=True)
    child = subprocess.Popen('wc -l %s' % vcffile, shell=True, stdout=subprocess.PIPE)
    total_line = int(child.communicate()[0].split()[0])
    print('total %s lines' % total_line)
    step = total_line / N
    print(1)
    cmd_first = "sed -n '1,%sp' %s > %s.1.vcf" % (step, vcffile, prefix)
    subprocess.call(cmd_first, shell=True)
    for i in range(2, N):
        print(i)
        st = (i - 1) * step + 1
        ed = i * step
        cmd = "sed -n '%s,%sp' %s > %s.%s.tmp.vcf" % (st, ed, vcffile, prefix, i)
        subprocess.call(cmd, shell=True)
    print(i + 1)
    cmd_last = "sed -n '%s,%sp' %s > %s.%s.tmp.vcf" % ((ed + 1), total_line, vcffile, prefix, (i + 1))
    subprocess.call(cmd_last, shell=True)
    for i in range(2, N + 1):
        cmd_cat = 'cat %s.header %s.%s.tmp.vcf > %s.%s.vcf' % (prefix, prefix, i, prefix, i)
        subprocess.call(cmd_cat, shell=True)

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
    
    

def merge_files(args):
    """
    %prog merge_files pattern out_fn
    combine split vcf files to a single one. Pattern example: 'hmp321_agpv4_chr9.%s.beagle.vcf'
    revise the lambda fucntion to fit your file patterns
    """

    p = OptionParser(merge_files.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    pattern,out_fn, = args

    fns = [str(i) for i in list(Path('.').glob(pattern))]
    fns_sorted = sorted(fns, key=lambda x: int(x.split('_')[1].split('-')[0]))
    print(fns_sorted)
    print('%s files found!'%len(fns_sorted))

    f = open(out_fn, 'w')
    print(fns_sorted[0])
    with open(fns_sorted[0]) as f1:
        for i in f1:
            f.write(i)
    for i in fns_sorted[1:]:
        print(i)
        with open(i) as f2:
            for j in f2:
                if not j.startswith('#'):
                    f.write(j)

def impute_beagle(args):
    """
    %prog impute_beagle dir_in dir_out
    impute missing data in vcf using beagle 
    """
    p = OptionParser(impute_beagle.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help = 'file pattern for vcf files in dir_in')
    p.add_option('--parameter_file',
                 help = 'file including window, overlap parameters')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')

    if opts.parameter_file:
        df = pd.read_csv(opts.parameter_file)
        df = df.set_index('chr')
    
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        print(vcf.name)
        chrom = int(vcf.name.split('.')[0].split('hr')[-1])
        print('chr : %s'%chrom)
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.BG'
        out_fn_path = out_path/out_fn
        if opts.parameter_file:
            window = df.loc[chrom, 'marker_10cm']
            overlap = df.loc[chrom, 'marker_2cm']
            print('window: %s; overlap: %s'%(window, overlap))
            cmd = 'beagle -Xmx60G gt=%s out=%s window=%s overlap=%s nthreads=10' % (vcf, out_fn_path, window, overlap)
        else:
            cmd = 'beagle -Xmx60G gt=%s out=%s nthreads=10' % (begle, vcf, out_fn_path)
        header = multiCPU_header % (10, 167, 65000, sm, sm, sm)
        header += 'ml beagle/4.1\n'
        header += cmd
        with open('%s.beagle.slurm' % sm, 'w') as f:
            f.write(header)

def vcf2hmp(args):
    """
    %prog vcf2hmp vcf
    convert vcf generated from beagle to hmp format using tassel
    """
    p = OptionParser(vcf2hmp.__doc__)
    p.set_slurm_opts(jn=True)
    p.add_option('--version', default='2', choices=('1', '2'),
                 help='specify the hmp type. 1: hyploid. 2: diploid')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    vcffile, = args
    prefix = '.'.join(vcffile.split('.')[0:-1])
    cmd = '%s -Xms512m -Xmx10G -fork1 -vcf %s -export -exportType HapmapDiploid\n' % (tassel, vcffile) \
        if opts.version == '2' \
        else '%s -Xms512m -Xmx10G -fork1 -vcf %s -export -exportType Hapmap\n' % (tassel, vcffile)
    header = Slurm_header % (opts.time, opts.memory, opts.prefix, opts.prefix, opts.prefix)
    header += 'module load java/1.8\n'
    header += cmd
    f = open('%s.vcf2hmp.slurm' % prefix, 'w')
    f.write(header)
    f.close()
    print('slurm file %s.vcf2hmp.slurm has been created, you can submit your job file.' % prefix)


def FixIndelHmp(args):
    """
    %prog FixIndelHmp hmp
    Fix the InDels problems in hmp file generated from Tassel.
    """
    p = OptionParser(FixIndelHmp.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    hmpfile, = args
    prefix = '.'.join(hmpfile.split('.')[0:-1])

    f = open(hmpfile)
    f1 = open('%s.Findel.hmp' % prefix, 'w')

    bict = {'A': 'T', 'T': 'C', 'C': 'A', 'G': 'A'}
    for i in f:
        j = i.split()
        if '+' in j[1]:
            nb = bict[j[1][0]] if j[1][0] != '+' else bict[j[1][-1]]
            tail = '\t'.join(j[11:]) + '\n'
            n_tail = tail.replace('+', nb)
            head = '\t'.join(j[0:11])
            n_head = head.replace('/+', '/%s' % nb) \
                if j[1][0] != '+' \
                else head.replace('+/', '%s/' % nb)
            n_line = n_head + '\t' + n_tail
            f1.write(n_line)
        elif '-' in j[1]:
            nb = bict[j[1][0]] if j[1][0] != '-' else bict[j[1][-1]]
            tail = '\t'.join(j[11:]) + '\n'
            n_tail = tail.replace('-', nb)
            head = '\t'.join(j[0:11])
            n_head = head.replace('/-', '/%s' % nb) \
                if j[1][0] != '-' \
                else head.replace('-/', '%s/' % nb)
            n_line = n_head + '\t' + n_tail
            f1.write(n_line)
        else:
            f1.write(i)
    f.close()
    f1.close()


def downsampling(args):
    pass


if __name__ == "__main__":
    main()
