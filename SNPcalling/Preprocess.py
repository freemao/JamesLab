# -*- coding: UTF-8 -*-

"""
Prepare fastq files ready for SNP calling
"""

import sys
import os.path as op
from pathlib import Path
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob
from schnablelab.apps.natsort import natsorted
import subprocess
from schnablelab.apps.headers import Slurm_header


def main():
    actions = (
        ('fastqc', 'check the reads quality'),
        ('trimmomatic', 'quality control on sample reads'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def fastqc(args):
    """
    %prog fastqc in_dir out_dir
        in_dir: the dir where fastq files are located
        out_dir: the dir saving fastqc reports

    generate slurm files for fastqc jobs
    """
    p = OptionParser(fastqc.__doc__)
    p.add_option("--pattern", default = '*.fastq', 
            help="the pattern of fastq files, qutation needed") 
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    fn_fq, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        print('%s does not exist...')
    dir_path = Path(in_dir)
    fqs = dir_path.glob(opts.pattern)
    for fq in fqs:
        prf = '.'.join(fq.name.split('.')[0:-1])
        print(prf)
        cmd = 'fastqc %s -o %s'%(str(fq), out_dir)
        header = Slurm_header%(10, 1, prf, prf, prf)
        header += 'ml fastqc\n'
        header += cmd
	with open('%s.fastqc.slurm'%(prf), 'w') as f:
            f.write(header)

def trimmomatic(args):
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


if __name__ == "__main__":
    main()
