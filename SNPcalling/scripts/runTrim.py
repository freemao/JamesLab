header = '''#!/bin/sh
#SBATCH --time=01:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=20000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load trimmomatic

'''

import os
from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('.fq.gz')]
#print 'Total %s fastq.gz files'%len(allfiles)
for i in allfiles:
    sm = i.split('.')[0]
    cmd1 = 'trimmomatic SE %s %s CROP:185 HEADCROP:5 SLIDINGWINDOW:4:15 MINLEN:30'%(i, sm+'.trimed.fq')
    cmd2 = 'gzip %s'%(sm+'.trimed.fq')
    jobfile = '%s.trimc.slurm'%sm
    f = open(jobfile, 'w')
    f.write(header%(sm,sm,sm))
    f.write(cmd1+'\n')
    f.write(cmd2+'\n')
    f.close()
    jobcmd = 'sbatch -p jclarke %s'%jobfile
    call(jobcmd, shell=True)
