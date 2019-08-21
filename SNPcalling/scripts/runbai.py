header = '''#!/bin/sh
#SBATCH --time=01:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=10000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load samtools/0.1

'''

import os
from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('sorted.bam')]
#print 'Total %s sorted.bam files'%len(allfiles)
for i in allfiles:
    SM = i.split('.')[0]
    cmd = 'samtools index %s'%i
    jobfile = '%s.index.slurm'%SM
    f = open(jobfile, 'w')
    f.write(header%(SM,SM,SM))
    f.write(cmd)
    f.close()
    jobcmd = 'sbatch -p jclarke %s'%jobfile
    #jobcmd = 'sbatch %s'%jobfile
    ##print cmd
    call(jobcmd, shell=True)
