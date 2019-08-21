header = '''#!/bin/sh
#SBATCH --time=50:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=80000       # Maximum memory required per CPU (in megabytes)
#SBATCH --job-name=%s
#SBATCH --error=./%s.err
#SBATCH --output=./%s.out

'''

import os
from subprocess import call

def main(ref, info):
    allfiles = [i for i in os.listdir('.') if i.endswith('sorted.bam')]
    #print 'Total %s sorted.bam files'%len(allfiles)
    f1 = open('bamfiles.fb.list', 'w')
    for i in allfiles :
        f1.write(i + '\n')
    f1.close()

    f2 = open(info)
    chrlist = [i.rstrip() for i in f2]
    for seq in chrlist:
        cmd = '/work/schnablelab/cmiao/SorghumGWAS/scripts/freebayes/bin/freebayes -r %s -f %s -C 1 -L bamfiles.fb.list > %s\n'%(seq, ref, "_".join(seq.split(':'))+'.vcf')
        #print cmd
        jobfile = '%s.fb.slurm'%("_".join(seq.split(':')))
        f = open(jobfile, 'w')
        f.write(header%(seq,seq,seq))
        f.write(cmd)
        f.close()
        jobcmd = 'sbatch -p jclarke %s'%jobfile
        #jobcmd = 'sbatch %s'%jobfile
        call(jobcmd, shell=True)

import sys
if len(sys.argv) == 3:
    main(*sys.argv[1:])
else:
    print('Path_ref ChrInfoFile')
