# -*- coding: UTF-8 -*-

"""
base class and functions for image processing
"""
import sys
import numpy as np
import pandas as pd
from PIL import Image
from tqdm import tqdm
from pathlib import Path
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

def main():
    actions = (
        ('CropBorder', 'crop borders'),
        ('BatchCropBorder', 'apply CropBorder on large number of images using HPC')
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def CropBorder(args):
    '''
    %prog crop_border img1 img2 img3 ...

    crop image borders using PIL. 
    If multiple images are provided, they will be cropped using the same crop size
    '''
    p = OptionParser(CropBorder.__doc__)
    p.add_option('--crop_dim', default = '410,0,2300,1675',
        help = 'the dimension (left,upper,right,lower) after cropping')
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    
    dim = ([int(i) for i in opts.crop_dim.split(',')])
    for img_fn in args:
        img_out_fn = Path(img_fn).name.replace('.png', '.crp.png')
        Image.open(img_fn).crop(dim).save(Path(opts.out_dir)/img_out_fn, 'PNG')
    
def BatchCropBorder(args):
    '''
    %prog in_dir out_dir

    apply CropBorder on a large number of images
    '''
    p = OptionParser(BatchCropBorder.__doc__)
    p.add_option('--pattern', default='*.png',
                 help="file pattern of png files under the 'dir_in'")
    p.add_option('--crop_dim', default = '410,0,2300,1675',
        help = 'the dimension (left,upper,right,lower) after cropping')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchCropBorder.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    pngs = in_dir_path.glob(opts.pattern)
    cmds = []
    for img_fn in pngs:
        cmd = 'python -m schnablelab.ImageProcessing.base CropBorder %s --crop_dim %s --out_dir %s'%(img_fn, opts.crop_dim, out_dir)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)
    
if __name__ == "__main__":
    main()