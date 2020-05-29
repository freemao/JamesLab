# -*- coding: UTF-8 -*-

"""
functions to process vcf files
"""
import sys
import pandas as pd
import os.path as op
import subprocess
from pathlib import Path
from subprocess import run
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, put2slurm

# the location of linkimpute, beagle executable
lkipt = op.abspath(op.dirname(__file__)) + '/../apps/LinkImpute.jar'
begle = op.abspath(op.dirname(__file__)) + '/../apps/beagle.24Aug19.3e8.jar'
tassel = op.abspath(op.dirname(__file__)) + '/../apps/tassel-5-standalone/run_pipeline.pl'

def main():
    actions = (
        ('BatchFilterMissing', 'apply FilterMissing on multiple vcf files'),
        ('BatchFilterMAF', 'apply FilterMissing on multiple vcf files'),
        ('BatchFilterHetero', 'apply FilterMissing on multiple vcf files'),
        ('IndexVCF', 'index vcf using bgzip and tabix'),
        ('splitVCF', 'split a vcf to several smaller files with equal size'),
        ('merge_files', 'combine split vcf or hmp files'),
        ('combineFQ', 'combine split fqs'),
        ('impute_beagle', 'impute vcf using beagle or linkimpute'),
        ('FixIndelHmp', 'fix the indels problems in hmp file converted from tassel'),
        ('FilterVCF', 'remove bad snps using bcftools'),
        ('only_ALT', 'filter number of ALT'),
        ('fixGTsep', 'fix the allele separator for beagle imputation'),
        ('SummarizeLD', 'summarize ld decay in log scale'),
        ('calculateLD', 'calculate r2 using Plink'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def BatchFilterMissing(args):
    """
    %prog in_dir

    apply FilterMissing on multiple vcf files
    """
    p = OptionParser(BatchFilterMissing.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help="file pattern of vcf files in the 'dir_in'")
    p.add_option('--missing_cutoff', default='0.7',
                 help='missing rate cutoff, SNPs higher than this cutoff will be removed')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchFilterMissing.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, = args
    in_dir_path= Path(in_dir)
    vcfs = in_dir_path.glob(opts.pattern)
    cmds = []
    for vcf in vcfs:
        cmd = "python -m schnablelab.SNPcalling.base FilterMissing %s --missing_cutoff %s"%(vcf, opts.missing_cutoff)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

def BatchFilterMAF(args):
    """
    %prog in_dir

    apply FilterMAF on multiple vcf files
    """
    p = OptionParser(BatchFilterMAF.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help="file pattern of vcf files in the 'dir_in'")
    p.add_option('--maf_cutoff', default='0.01',
                 help='maf cutoff, SNPs lower than this cutoff will be removed')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchFilterMAF.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, = args
    in_dir_path= Path(in_dir)
    vcfs = in_dir_path.glob(opts.pattern)
    cmds = []
    for vcf in vcfs:
        cmd = "python -m schnablelab.SNPcalling.base FilterMAF %s --maf_cutoff %s"%(vcf, opts.maf_cutoff)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

def BatchFilterHetero(args):
    """
    %prog in_dir

    apply FilterMAF on multiple vcf files
    """
    p = OptionParser(BatchFilterHetero.__doc__)
    p.add_option('--pattern', default='*.vcf',
                 help="file pattern of vcf files in the 'dir_in'")
    p.add_option('--het_cutoff', default='0.1',
                 help='heterozygous rate cutoff, SNPs higher than this cutoff will be removed')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchFilterHetero.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, = args
    in_dir_path= Path(in_dir)
    vcfs = in_dir_path.glob(opts.pattern)
    cmds = []
    for vcf in vcfs:
        cmd = "python -m schnablelab.SNPcalling.base FilterHetero %s --het_cutoff %s"%(vcf, opts.het_cutoff)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

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

def calculateLD(args):
    """
    %prog vcf_fn/plink_prefix genome_size(bp) num_SNPs

    calculate LD using Plink
    args:
        vcf_fn/plink_prefix: specify either vcf/vcf.gz file or the prefix of plink bed/bim/fam files. 
        genome_size(bp): the size of the reference genome in bp. For reference: sorghum 68,400,000
        num_SNPs: the number of SNPs in the genotype file.
    """
    p = OptionParser(calculateLD.__doc__)
    p.add_option('--maf_cutoff', default='0.01',
                 help='only use SNP with the MAF higher than this cutoff to calculate LD')
    p.add_option('--max_distance', type='int', default=1000000,
                 help='the maximum distance for a pair of SNPs to calcualte LD')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=calculateLD.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_fn, g_size, n_snps, = args
    in_fn, g_size, n_snps = Path(in_fn), int(g_size), int(n_snps)

    if in_fn.name.endswith('.vcf') or in_fn.name.endswith('.vcf.gz'):
        input = f'--vcf {in_fn}'
    else:
        input = f'--bfile {in_fn}'
    n = 10
    ld_window, ld_window_bp = [], [] 
    while True:
        ld_window.append(n)
        dist = g_size//n_snps*n
        ld_window_bp.append(dist)
        n *= 10
        if dist>=1000000:
            break
    
    cmds = []
    for win_snp, win_bp in zip(ld_window, ld_window_bp):
        out_fn = Path(in_fn).name.split('.')[0]
        prob = 10/win_snp
        cmd = f'plink {input} --thin {prob} --r2 --ld-window 10 --ld-window-kb {win_bp//1000} --ld-window-r2 0 --maf {opts.maf_cutoff} --out {out_fn}.r2.thin{prob}.csv'
        cmds.append(cmd)
        print(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml plink'
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['cmd_header'] = cmd_header
        put2slurm(cmds, put2slurm_dict)

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
