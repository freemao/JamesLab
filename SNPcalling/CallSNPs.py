# -*- coding: UTF-8 -*-

"""
Call SNPs on HTS data using GATK, Freebayes.
"""

import os
import sys
import numpy as np
import pandas as pd
import os.path as op
from pathlib import Path
from subprocess import run
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm


def main():
    actions = (
        ('genGVCF', 'generate gvcf for each sample using GATK HaplotypeCaller'),
        ('freebayes', 'call SNPs using freebayes'),
        ('samtools', 'call SNPs using samtools'),
        ('gatk', 'call SNPs using gatk'),
)
    p = ActionDispatcher(actions)
    p.dispatch(globals())
    
def genGVCF(args):
    """
    %prog genGVCF ref.fa bams.csv region.txt out_dir

    run GATK HaplotypeCaller in GVCF mode
    args:
        ref.fa: reference sequence file
        bams.csv: csv file containing all bam files and their sample names
        region.txt: genomic intervals defined by each row to speed up GVCF calling. 
            example regions: Chr01, Chr01:1-100
        out_dir: where the gVCF files save to
    """
    p = OptionParser(genGVCF.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    ref, bams_csv, region_txt, out_dir, = args
    out_dir_path = Path(out_dir)
    if not out_dir_path.exists():
        sys.exit(f'output directory {out_dir_path} does not exist!')
    
    regions = []
    with open(region_txt) as f:
        for i in f:
            regions.append(i.rstrip())
    
    mem = int(opts.memory)//1024

    print('defined genomic intervals: %s'%(','.join(regions)))
    df_bam = pd.read_csv(bams_csv)
    cmds = []
    for sm, grp in df_bam.groupby('sm'):
        print(f'{grp.shape[0]} bam files for sample {sm}')
        input_bam = '-I ' + ' -I '.join(grp['fn_path'].tolist())
        output_fn = f'{sm}.g.vcf'
        for region in regions:
            print(f'region: {region}')
            # --sample-name: Name of single sample to use from a multi-sample bam
            cmd = f"gatk --java-options '-Xmx{mem}g' HaplotypeCaller -R {ref} {input_bam} -O {out_dir_path/output_fn} --sample-name {sm} --emit-ref-confidence GVCF -L {region.rstrip()}"
            cmds.append(cmd)
    
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print(f'check {cmd_sh} for all the commands!')

    cmd_header = 'ml gatk4'
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm_dict['cmd_header'] = cmd_header
        put2slurm(cmds, put2slurm_dict)

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
            if ':0-' in reg:
                reg = reg.replace(':0-', ':1-')
            reg_fn = reg.replace(':','_')
            reg_fn_vcf = '%s.gatk.vcf'%reg_fn
            reg_fn_vcf_path = out_path/reg_fn_vcf
            cmd = "gatk --java-options '-Xmx13G' HaplotypeCaller \\\n-R %s -L %s \\\n%s-O %s"%(ref, reg, inputs, reg_fn_vcf_path)
            header = Slurm_header%(165, 15000, reg_fn, reg_fn, reg_fn)
            header += 'ml gatk4/4.1\n'
            header += cmd
            with open('%s.gatk.slurm'%reg_fn, 'w') as f1:
                f1.write(header)

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
            cmd = 'freebayes -r %s -f %s -C 1 -F 0.05 -L %s -u -n 2 -g %s > %s\n'%(reg, ref, bams,opts.max_depth, reg_fn_vcf_path)
            header = Slurm_header%(165, 50000, reg_fn, reg_fn, reg_fn)
            header += 'ml freebayes/1.3\n'
            header += cmd
            with open('%s.fb.slurm'%reg_fn, 'w') as f1:
                f1.write(header)
            print('slurm files %s.fb.slurm has been created'%reg_fn)

if __name__ == "__main__":
    main()
