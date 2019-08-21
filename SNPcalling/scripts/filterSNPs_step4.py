header = '''#!/bin/sh
#SBATCH --time=01:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load bcftools

'''

import os
from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('.step3.vcf')]
#allfiles = [i for i in os.listdir('.') if i.endswith('.step1.vcf')]
print('Total %s files'%len(allfiles))
for i in allfiles:
    prefix = i.split('.')[0]
    new_fn = prefix + '.step4.vcf' 
    cmd = "python /work/schnablelab/cmiao/SorghumGWAS/scripts/Rm_BadIndels.py %s %s"%(i, new_fn)
    jobfile = '%s.step4.slurm'%prefix
    f = open(jobfile, 'w')
    f.write(header%(prefix, prefix, prefix))
    f.write(cmd)
    f.close()
    jobcmd = 'sbatch -p jclarke %s'%jobfile
    #jobcmd = 'sbatch %s'%jobfile
    call(jobcmd, shell=True)
