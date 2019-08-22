from schnablelab.apps.base import glob
from pathlib import Path
import sys
from subprocess import run

pfx_r1 = glob('*_trim.R1.fastq')
pfx_r2 = glob('*_trim.R2.fastq')
pfx_un = glob('*_trim.Unpaired.fastq')

for r1,r2,un in zip(pfx_r1, pfx_r2, pfx_un):
    sm1 = r1.split('_trim')[0]
    sm2 = r2.split('_trim')[0]
    if sm1 == sm2:
        cmd = 'python -m schnablelab.SNPcalling.Preprocess align /work/schnablelab/cmiao/TimeSeriesGWAS/Genotype_GBS/Reference_Genome_4th/Sbicolor_454_v3.0.1 %s %s'%(r1, r2)
        run(cmd, shell=True)
    else:
        sys.exit('sm1 != sm2')
    cmd1 = 'python -m schnablelab.SNPcalling.Preprocess align /work/schnablelab/cmiao/TimeSeriesGWAS/Genotype_GBS/Reference_Genome_4th/Sbicolor_454_v3.0.1 %s'%un
    run(cmd1, shell=True)
