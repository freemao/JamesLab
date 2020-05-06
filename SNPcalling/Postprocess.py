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
import numpy as np
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
        ('IndexVCF', 'index vcf using bgzip and tabix'),
        ('splitVCF', 'split a vcf to several smaller files with equal size'),
        ('merge_files', 'combine split vcf or hmp files'),
        ('combineFQ', 'combine split fqs'),
        ('impute_beagle', 'impute vcf using beagle or linkimpute'),
        ('vcf2hmp', 'convert vcf to hmp format'),
        ('FixIndelHmp', 'fix the indels problems in hmp file converted from tassel'),
        ('FilterVCF', 'remove bad snps using bcftools'),
        ('only_Missing', 'perform filtering only on missing data'),
        ('only_ALT', 'filter number of ALT'),
        ('only_MAF', 'filter MAF'),
        ('only_Hetero', 'filter high heterozygous loci'),
        ('fixGTsep', 'fix the allele separator for beagle imputation'),
        ('SummarizeLD', 'summarize ld decay in log scale'),
        ('EstimateLD', 'estimate ld using tassel'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def only_MAF(args):
    """
    %prog in_dir out_dir

    filter MAF
    """
    p = OptionParser(only_MAF.__doc__)
    p.set_slurm_opts(jn=True)
    p.add_option('--pattern', default='*.vcf',
                 help='file pattern for vcf files in dir_in')
    p.add_option('--maf', default='0.01',
                 help='maf cutoff')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcffile in vcfs:
        prefix = '.'.join(vcffile.name.split('.')[0:-1])
        cmd = "python -m schnablelab.SNPcalling.base MAF %s %s\n"%(vcffile, opts.maf)
        with open('%s.maf.slurm'%prefix, 'w') as f:
            header = Slurm_header%(opts.time, opts.memory, prefix, prefix, prefix)
            header += 'ml bcftools\n'
            header += cmd
            f.write(header)
            print('slurm file %s.maf.slurm has been created, you can sbatch your job file.'%prefix)

def only_ALT(args):
    """
    %prog in_dir out_dir

    filter number of ALT using bcftools
    """
    p = OptionParser(only_ALT.__doc__)
    p.set_slurm_opts(jn=True)
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
    for vcffile in vcfs:
        prefix = '.'.join(vcf.name.split('.')[0:-1])
        new_f = prefix + '.alt1.vcf'
        cmd = "bcftools view -i 'N_ALT=1' %s > %s"%(vcffile, new_f)
        with open('%s.alt1.slurm'%prefix, 'w') as f:
            header = Slurm_header%(opts.time, opts.memory, prefix, prefix, prefix)
            header += 'ml bacftools\n'
            header += cmd
            f.write(header)
            print('slurm file %s.alt1.slurm has been created, you can sbatch your job file.'%prefix)

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

def IndexVCF(args):
    """
    %prog IndexVCF in_dir out_dir

    index vcf using bgzip and tabix
    """
    p = OptionParser(IndexVCF.__doc__)
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
        out_fn = vcf.name+'.gz'
        out_fn_path = out_path/out_fn
        cmd1 = 'bgzip -c %s > %s\n'%(vcf, out_fn_path)
        cmd2 = 'tabix -p vcf %s\n'%(out_fn_path)
        header = Slurm_header%(10, 20000, sm, sm, sm)
        header += 'ml tabix\n'
        header += cmd1
        header += cmd2
        with open('%s.idxvcf.slurm'%sm, 'w') as f:
            f.write(header)

def only_Hetero(args):
    """
    %prog only_Hetero in_dir out_dir

    filter SNPs with high heterozygous rate
    """

    p = OptionParser(only_Hetero.__doc__)
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
        cmd = 'python -m schnablelab.SNPcalling.base Heterozygous %s %s --h2_rate %s'%(vcf, out_path/out_fn, opts.rate)
        header = Slurm_header%(10, 20000, sm, sm, sm)
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
    subprocess.call('module load bcftools', shell=True)
    child = subprocess.Popen('bcftools query -l %s|wc -l'%vcffile, shell=True, stdout=subprocess.PIPE)
    SMs_num = int(child.communicate()[0])
    return SMs_num


def only_Missing(args):
    """
    %prog dir_in
    Remove SNPs with high missing rate (>0.7 by default)
    """
    p = OptionParser(only_Missing.__doc__)
    p.add_option('--pattern', default = '*.vcf',
        help = 'remove loci with high missing data rate')
    p.add_option('--missing_rate', default = 0.7,
        help = 'specify the missing rate cutoff. SNPs with missing rate higher than this cutoff will be removed.')

    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, = args
    dir_path = Path(in_dir)
    vcfs = dir_path.glob(opts.pattern)
    for vcf in vcfs:
        prefix = '.'.join(vcf.name.split('.')[0:-1])
        cmd = "python -m schnablelab.SNPcalling.base Missing %s --missing_rate %s\n"%(vcf, float(opts.missing_rate))
        header = Slurm_header%(10, 5000, prefix, prefix, prefix)
        header += 'ml bcftools\n'
        header += cmd
        with open('%s.mis%s.slurm'%(prefix, opts.missing_rate), 'w') as f:
            f.write(header)
        print('%s.mis%s.slurm has been generated!'%(prefix, opts.missing_rate))

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
    fns_sorted = sorted(fns, key=lambda x: int(x.split('.')[0][3:]))
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
    cmd = 'run_pipeline.pl -Xms512m -Xmx10G -fork1 -vcf %s -export -exportType HapmapDiploid\n'%vcffile \
        if opts.version == '2' \
        else 'run_pipeline.pl -Xms512m -Xmx10G -fork1 -vcf %s -export -exportType Hapmap\n'%vcffile
    header = Slurm_header % (opts.time, opts.memory, opts.prefix, opts.prefix, opts.prefix)
    header += 'module load java/1.8\n'
    header += 'module load tassel/5.2\n'
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

def EstimateLD(args):
    """
    %prog dir_in dir_out
    run LD decay using tassel
    """
    p = OptionParser(EstimateLD.__doc__)
    p.set_slurm_opts(jn=True)
    p.add_option('--pattern', default='*vcf',
                 help='pattern of vcf files')
    p.add_option('--window_size', default='1000',
                 help='specify how many SNPs in the sliding window')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    dir_in, dir_out = args
    dir_out = Path(dir_out)
    if not dir_out.exists():
        dir_out.mkdir()
    for vcf in Path(dir_in).glob(opts.pattern):
        prefix = vcf.name.replace('.vcf', '')
        out_fn = '%s.ld'%prefix
        cmd = 'run_pipeline.pl -Xms512m -Xmx14g -fork1 -vcf %s -ld -ldWinSize %s -ldType SlidingWindow -td_tab %s/%s\n'%(vcf, opts.window_size, dir_out, out_fn)
        header = Slurm_header % (opts.time, 15000, prefix, prefix, prefix)
        header += 'ml java/1.8\n'
        header += 'ml tassel/5.2\n'
        header += cmd
        with open('%s.estLD.slurm' % prefix, 'w') as f:
            f.write(header)
        print('slurm file %s.estLD.slurm has been created, you can submit your job file.' % prefix)

def SummarizeLD(args):
    """
    %prog dir_in dir_out
    summarize LD decay in log scale
    """
    p = OptionParser(EstimateLD.__doc__)
    p.set_slurm_opts(jn=True)
    p.add_option('--pattern', default='*.ld.txt',
                 help='pattern of ld.txt files')
    p.add_option('--max_dist', default='1,000,000',
                 help='the maximum ld distance')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    dir_in, dir_out = args
    dir_out = Path(dir_out)
    if not dir_out.exists():
        dir_out.mkdir()
    num0 = opts.max_dist.count('0')

    for fn in Path(dir_in).glob(opts.pattern):
        prefix = '.'.join(fn.name.split('.')[0:-1])
        out_fn = '%s.sum.csv'%prefix
        cmd = 'python -m schnablelab.SNPcalling.base SummarizeLD %s %s %s/%s\n'%(fn, num0, dir_out, out_fn)
        header = Slurm_header % (opts.time, opts.memory, prefix, prefix, prefix)
        header += cmd
        with open('%s.sumLD.slurm' % prefix, 'w') as f:
            f.write(header)
        print('slurm file %s.sumLD.slurm has been created, you can submit your job file.' % prefix)


if __name__ == "__main__":
    main()
