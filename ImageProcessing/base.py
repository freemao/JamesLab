# -*- coding: UTF-8 -*-

"""
base class and functions for image processing
"""
import csv
import cv2
import sys
import numpy as np
import pandas as pd
from PIL import Image
from tqdm import tqdm
from pathlib import Path
from skimage.util import invert
from skimage.morphology import convex_hull_image
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

def main():
    actions = (
        ('CropFrame', 'crop borders'),
        ('BatchCropFrame', 'apply CropBorder on large number of images on HPC'),
        ('CropObject', 'crop based on object box'),
        ('BatchCropObject', 'apply CropObject on large number of images on HPC'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

class ProsImage():
    '''
    basic functions of processing image
    '''
    def __init__(self, filename):
        self.fn = filename
        # four channels: RGBA
        self.PIL_img = Image.open(filename).convert('RGB')
        self.width, self.height = self.PIL_img.size
        self.array = np.array(self.PIL_img)
        self.array_g = self.array[:,:,1]
        # 2*green/(red+blue)
        #self.array_gidx = (2*self.array_g)/(self.array[:,:,0]+self.array[:,:,2])
        self.h_w_c = self.array.shape

    def crop(self, crp_dim):
        return self.PIL_img.crop(crp_dim)
    
    def binary_thresh(self, method='green_channel'):
        if method == 'green_channel':
            # background: white(255), plant: black(0)
            _, thresh = cv2.threshold(self.array_g, 130, 255, cv2.THRESH_BINARY)
        elif method == 'green_index':
            _, thresh = cv2.threshold(self.array_gidx, 1.12, 255, cv2.THRESH_BINARY)
        # background: black(0), plant: whilte(255)
        thresh_ivt = invert(thresh)
        return thresh_ivt

    def box_size(self, thresh_method='green_channel'):
        thresh_ivt = self.binary_thresh(method=thresh_method)
        height = np.sum(thresh_ivt, axis=1)
        width = np.sum(thresh_ivt, axis=0)
        idx_top = next(x for x, val in enumerate(height) if val > 0)
        idx_bottom = self.height-next(x for x, val in enumerate(height[::-1]) if val > 0)
        idx_left = next(x for x, val in enumerate(width) if val > 0) 
        idx_right = self.width-next(x for x, val in enumerate(width[::-1]) if val > 0) 
        return (idx_top, idx_bottom, idx_left, idx_right)

    def hull(self, thresh_method='green_channel'):
        thresh_ivt = self.binary_thresh(method=thresh_method)        
        # non-hull part: False, hull part: True
        chull = convex_hull_image(thresh_ivt)
        # view the hull image
        #plant part: 2, non-hull part: 0 (False) (black), hull-part excluding plant: 1 (True) gray
        chull_diff = np.where(thresh_ivt == 255, 2, chull)

def CropObject(args):
    '''
    %prog CropObject img1 img2 img3 ...

    crop image based the threshold box
    '''
    p = OptionParser(CropObject.__doc__)
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    p.add_option('--pad', type='int', default=50,
        help = 'specify the pad size')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    
    for img_fn in args:
        img_out_fn = Path(img_fn).name.replace('.png', '.CrpObj.png')
        img = ProsImage(img_fn)

        idx_top, idx_bottom, idx_left, idx_right = img.box_size()
        print('box size:\n  top:%s, bottom:%s, left:%s, right:%s'%(idx_top, idx_bottom, idx_left, idx_right))
        
        idx_top -= opts.pad
        idx_bottom += opts.pad
        idx_left -= opts.pad
        idx_right += opts.pad

        idx_top = 0 if idx_top < 0 else idx_top
        idx_bottom = img.height if idx_bottom > img.height else idx_bottom
        idx_left = 0 if idx_left < 0 else idx_left
        idx_right = img.width if idx_right > img.width else idx_right
        crp_dim = (idx_left, idx_top, idx_right, idx_bottom)
        img.crop(crp_dim).save(Path(opts.out_dir)/img_out_fn, 'PNG')

def BatchCropObject(args):
    '''
    %prog in_dir out_dir

    apply BatchCropObject on a large number of images
    '''
    p = OptionParser(BatchCropObject.__doc__)
    p.add_option('--pattern', default='*.png',
                 help="file pattern of png files under the 'dir_in'")
    p.add_option('--pad', type='int', default=50,
        help = 'specify the pad size')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchCropObject.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    pngs = in_dir_path.glob(opts.pattern)
    cmds = []
    for img_fn in pngs:
        cmd = 'python -m schnablelab.ImageProcessing.base CropObject %s --out_dir %s --pad %s'%(img_fn, out_dir, opts.pad)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)


def CropFrame(args):
    '''
    %prog CropFrame img1 img2 img3 ...

    crop image borders using PIL. 
    If multiple images are provided, they will be cropped using the same crop size
    '''
    p = OptionParser(CropFrame.__doc__)
    p.add_option('--crop_dim', default = '410,0,2300,1675',
        help = 'the dimension (left,upper,right,lower) after cropping')
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    
    dim = ([int(i) for i in opts.crop_dim.split(',')])
    for img_fn in args:
        img_out_fn = Path(img_fn).name.replace('.png', '.CrpFrm.png')
        img = ProsImage(img_fn)
        img.crop(dim).save(Path(opts.out_dir)/img_out_fn, 'PNG')
    
def BatchCropFrame(args):
    '''
    %prog in_dir out_dir

    apply BatchCropFrame on a large number of images
    '''
    p = OptionParser(BatchCropFrame.__doc__)
    p.add_option('--pattern', default='*.png',
                 help="file pattern of png files under the 'dir_in'")
    p.add_option('--crop_dim', default = '410,0,2300,1675',
        help = 'the dimension (left,upper,right,lower) after cropping')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchCropFrame.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    pngs = in_dir_path.glob(opts.pattern)
    cmds = []
    for img_fn in pngs:
        cmd = 'python -m schnablelab.ImageProcessing.base CropFrame %s --crop_dim %s --out_dir %s'%(img_fn, opts.crop_dim, out_dir)
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)
    
if __name__ == "__main__":
    main()