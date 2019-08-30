# -*- coding: UTF-8 -*-

"""
Split a big file using sed.
Find more details at cnblog:
www.cnblogs.com/freemao/p/7076127.html
"""

from pathlib import Path
import os.path as op
import sys
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob
from schnablelab.apps.natsort import natsorted
import subprocess
from subprocess import run
from schnablelab.apps.headers import Slurm_header

# the location of linkimpute, beagle executable
lkipt = op.abspath(op.dirname(__file__)) + '/../apps/LinkImpute.jar'
begle = op.abspath(op.dirname(__file__)) + '/../apps/beagle.24Aug19.3e8.jar'
tassel = op.abspath(op.dirname(__file__)) + '/../apps/tassel-5-standalone/run_pipeline.pl'


def main():
    actions = (
        ('splitVCF', 'split a vcf to several smaller files with equal size'),
        ('combineVCF', 'combine split vcfs'),
        ('combineFQ', 'combine split fqs'),
        ('impute', 'impute vcf using beagle or linkimpute'),
        ('vcf2hmp', 'convert vcf to hmp format'),
        ('FixIndelHmp', 'fix the indels problems in hmp file converted from tassel'),
        ('FilterVCF', 'remove bad snps using bcftools'),
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
        out_fn_path = out_path/out_fn
        cmd = 'python -m schnablelab.SNPcalling.FilterSNPs Heterozygous %s %s --h2_rate %s'%(vcf, out_fn_path, opts.rate)
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
    p.add_option('--qual', default='10',
                 help='minimum snp quality, 10: 10% is wrong, 20: 1% is wrong')
    p.add_option('--n_alt', default='1',
                 help='number of alt')
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

    cond1 = 'N_ALT==%s && QUAL>=%s'%(opts.n_alt, opts.qual)
    if opts.maf:
        cond1 += ' && MAF>=%s'%opts.maf
    if opts.missing:
        cond1 += ' && NS/N_SAMPLES > %s'%opts.missing
    cmd = "bcftools view -i '{cond1}' -v '{stype}' %s".format(cond1=cond1, stype=opts.stype)
    if opts.normalization:
        cmd += ' | bcftools norm -f %s -m -both'%(opts.ref)
    cmd += ' > %s'
    
    for vcf in vcfs:
        sm = '.'.join(vcf.name.split('.')[0:-1])
        out_fn = sm+'.bcflt.vcf'
        out_fn_path = out_path/out_fn
        header = Slurm_header%(10, 20000, sm, sm, sm)
        header += 'ml bcftools\n'
        header += cmd%(vcf, out_fn_path)
        with open('%s.bcflt.slurm'%sm, 'w') as f:
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
    
    

def combineVCF(args):
    """
    %prog combineVCF N pattern
    combine split vcf (1-based) files to a single one. Pattern example: hmp321_agpv4_chr9.%s.beagle.vcf
    """

    p = OptionParser(combineVCF.__doc__)
    p.add_option('--header', default='yes', choices=('yes', 'no'),
                 help='choose whether add header or not')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    N, vcf_pattern, = args
    N = int(N)
    new_f = vcf_pattern.replace('%s', '').replace('..', '.')
    print('output file: %s' % new_f)

    f = open(new_f, 'w')

    fn1 = open(vcf_pattern % 1)
    print(1)
    if opts.header == 'yes':
        for i in fn1:
            f.write(i)
    else:
        for i in fn1:
            if not i.startswith('#'):
                f.write(i)
    fn1.close()
    for i in range(2, N + 1):
        print(i)
        fn = open(vcf_pattern % i)
        for j in fn:
            if not j.startswith('#'):
                f.write(j)
        fn.close()
    f.close()


def impute(args):
    """
    %prog impute dir_in dir_out
    impute missing data in vcf using beagle or linkimpute
    """
    p = OptionParser(impute.__doc__)
    p.add_option('--software', default='linkimpute', choices=('linkimpute', 'beagle'),
                 help='specify the imputation software')
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
        out_fn = sm+'.LK.vcf' if opts.software=='linkimpute' else sm+'.BG'
        out_fn_path = out_path/out_fn
        cmd = 'java -Xmx60G -jar %s -v %s %s' % (lkipt, vcf, out_fn_path) \
            if opts.software == 'linkimpute' \
            else 'java -Xmx60G -jar %s gt=%s out=%s' % (begle, vcf, out_fn_path)
        header = Slurm_header % (165, 61000, sm, sm, sm)
        header += 'ml java/1.7\n' if opts.software == 'linkimpute' else 'ml java/1.8\n'
        header += cmd
        with open('%s.%s.slurm' % (sm, opts.software), 'w') as f:
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
