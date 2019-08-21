header = '''#!/bin/sh
#SBATCH --time=05:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=20000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

module load bwa/0.7

'''

import os
import os.path

from subprocess import call
allfiles = [i for i in os.listdir('.') if i.endswith('trimed.fq.gz') and not os.path.isfile(i.split('.')[0]+'.sam')]
#print allfiles
#print 'Total %s fq.gz files'%len(allfiles)
def runbwa(refversion):
    for i in allfiles:
        SM = i.split('.')[0]
        R = r"'@RG\tID:%s\tSM:%s'"%(SM, SM)
        samoutput = '%s.sam'%SM
        if refversion == '3':
            cmd = 'bwa mem -R %s %s %s > %s \n'%(R, '/work/schnablelab/cmiao/SorghumGWAS/shared_files/references/Sbicolor_313_v3.0.fa', i, samoutput)
        elif refversion == '1':
            cmd = 'bwa mem -R %s %s %s > %s \n'%(R, '/work/schnablelab/cmiao/SorghumGWAS/shared_files/references/Sbicolor_79.fa', i, samoutput)
        else:
            print('only 1 or 3')
        jobfile = '%s.bwa.slurm'%SM
        f = open(jobfile, 'w')
        f.write(header%(SM,SM,SM))
        f.write(cmd)
        f.close()
        jobcmd = 'sbatch -p jclarke %s'%jobfile
        call(jobcmd, shell=True)
import sys
if len(sys.argv)==2:
    runbwa(*sys.argv[1:])
else:
    print('python runbwa.py 1 or 3')
