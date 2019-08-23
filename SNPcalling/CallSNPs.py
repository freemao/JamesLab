# -*- coding: UTF-8 -*-

"""
Call SNPs on GBS data using freebayes.
"""

import os
import os.path as op
import sys
import pandas as pd
import numpy as np
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
    %prog freebayes region.txt ref.fa bam_list.txt 

    create freebayes slurm jobs for each splitted region defined in region.txt file
    """
    p = OptionParser(freebayes.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    region, ref, bams, = args

    with open(region) as f:
        for reg in f:
            reg = reg.strip()
            reg_fn = reg.replace(':','_')
            reg_fn_vcf = '%s.vcf'%reg_fn
            cmd = 'freebayes -r %s -f %s -C 1 -F 0.05 -L %s -u -n 2 > %s\n'%(reg, ref, bams, reg_fn_vcf)
            header = Slurm_header%(165, 20000, reg_fn, reg_fn, reg_fn)
            header += 'ml freebayes\n'
            header += cmd
            with open('%s.fb.slurm'%reg_fn, 'w') as f1:
                f1.write(header)
            print('slurm files %s.fb.slurm has been created'%reg_fn)

if __name__ == "__main__":
    main()
