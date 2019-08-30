# -*- coding: UTF-8 -*-

"""
Call SNPs on GBS data using freebayes.
"""

import os
import os.path as op
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from schnablelab.apps.base import ActionDispatcher, OptionParser
from schnablelab.apps.headers import Slurm_header
from subprocess import run

def main():
    actions = (
        ('freebayes', 'call SNPs using freebayes'),
        ('samtools', 'call SNPs using samtools'),
        ('gatk', 'call SNPs using gatk'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())
    
def freebayes(args):
    """
    %prog freebayes region.txt ref.fa bam_list.txt out_dir

    create freebayes slurm jobs for each splitted region defined in region.txt file
    """
    p = OptionParser(freebayes.__doc__)
    p.add_option('--max_depth', default=10000,
            help = 'cites where the mapping depth higher than this value will be ignored')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    region, ref, bams,out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')

    with open(region) as f:
        for reg in f:
            reg = reg.strip()
            reg_fn = reg.replace(':','_')
            reg_fn_vcf = '%s.fb.vcf'%reg_fn
            reg_fn_vcf_path = out_path/reg_fn_vcf
            cmd = 'freebayes -r %s -f %s -C 1 -F 0.05 -L %s -u -n 1 -g %s > %s\n'%(reg, ref, bams,opts.max_depth, reg_fn_vcf_pth)
            header = Slurm_header%(165, 50000, reg_fn, reg_fn, reg_fn)
            header += 'ml freebayes/1.3\n'
            header += cmd
            with open('%s.fb.slurm'%reg_fn, 'w') as f1:
                f1.write(header)
            print('slurm files %s.fb.slurm has been created'%reg_fn)

def gatk(args):
    """
    %prog gatk ref.fa bam_list.txt region.txt out_dir

    run GATK HaplotypeCaller
    """
    p = OptionParser(gatk.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref, bams, regions, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    with open(bams) as f:
        inputs = ''.join(['-I %s \\\n'%(i.rstrip()) for i in f])
    with open(regions) as f:
        for reg in f:
            reg = reg.strip()
            reg_fn = reg.replace(':','_')
            reg_fn_vcf = '%s.gatk.vcf'%reg_fn
            reg_fn_vcf_path = out_path/reg_fn_vcf
            cmd = "gatk --java-options '-Xmx13G' HaplotypeCaller \\\n-R %s \\\n%s-O %s"%(ref, inputs, reg_fn_vcf_path)
            header = Slurm_header%(165, 15000, reg_fn, reg_fn, reg_fn)
            header += 'ml gatk4/4.1\n'
            header += cmd
            with open('%s.gatk.slurm'%reg_fn, 'w') as f1:
                f1.write(header)

if __name__ == "__main__":
    main()
