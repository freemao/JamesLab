header = '''#!/bin/sh
#SBATCH --time=01:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=20000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load samtools/0.1

'''

import os
from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('bam')]
#print 'Total %s bam files'%len(allfiles)
for i in allfiles:
    SM = i.split('.')[0]
    output = '%s.sorted'%SM
    cmd = 'samtools sort %s %s'%(i, output)
    jobfile = '%s.sortbam.slurm'%SM
    f = open(jobfile, 'w')
    f.write(header%(SM,SM,SM))
    f.write(cmd)
    f.close()
    jobcmd = 'sbatch -p jclarke %s'%jobfile
    #jobcmd = 'sbatch  %s'%jobfile
    call(jobcmd, shell=True)
