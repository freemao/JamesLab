header = '''#!/bin/sh
#SBATCH --time=05:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=10000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load bcftools

'''

import os
import os.path

from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('.vcf')]
print(allfiles)
print('Total %s .vcf files'%len(allfiles))
for i in allfiles:
    SM = i.split('.')[0]
    cmd = "bcftools view -i 'N_ALT==1 && QUAL>=10 && MAF>=0.01 && NS/N_SAMPLES > 0.2' -v 'snps,indels' %s | bcftools norm -f /work/schnablelab/cmiao/SorghumGWAS/shared_files/references/Sbicolor_79.fa -m -both > %s.flt.vcf"%(i, SM)
    print(cmd)
    jobfile = '%s.bcftools.slurm'%SM
    f = open(jobfile, 'w')
    f.write(header%(SM,SM,SM))
    f.write(cmd)
    f.close()
    jobcmd = 'sbatch -p jclarke %s'%jobfile
    call(jobcmd, shell=True)
