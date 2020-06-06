# -*- coding: UTF-8 -*-

"""
base class and functions for image processing
"""
import csv
import cv2
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image
from tqdm import tqdm
from pathlib import Path
from skimage import io
from skimage.util import invert
from collections import namedtuple
from skimage.morphology import convex_hull_image
from matplotlib.patches import Arrow, Circle
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

def main():
    actions = (
        ('PlantArea', 'calculate plant pixel area'),
        ('toJPG', 'convert png to jpg'),
        ('Batch2JPG', 'apply toJPG on HPC'),
        ('Resize', 'resize images'),
        ('BatchResize', 'apply Resize on large number of images on HPC'),
        ('CropFrame', 'crop borders'),
        ('BatchCropFrame', 'apply CropBorder on large number of images on HPC'),
        ('CropObject', 'crop based on object box'),
        ('BatchCropObject', 'apply CropObject on large number of images on HPC'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def ShowMarkers(img_fn, coordinate_array, colors):
    '''
    show markers on an image
    args:
        img_fn: image filename
        coordinate_array: 2D array with shape (, 2)
        colors: list of colors for each coordinate
    return:
        matplotlib ax with patches added 
    '''
    img = io.imread(img_fn)
    _, ax = plt.subplots(1)
    ax.imshow(img)
    for xy, c in zip(coordinate_array, colors):
        circle = Circle(xy, radium=2, color=c)
        ax.add_patch(circle)
    return ax

class ProsImage():
    '''
    basic functions of processing image
    '''
    def __init__(self, filename):
        self.fn = filename
        self.format = filename.split('.')[-1]
        self.PIL_img = Image.open(filename)
        self.width, self.height = self.PIL_img.size
        self.array = np.array(self.PIL_img)[:,:,0:3]
        self.array_r = self.array[:,:,0]
        self.array_g = self.array[:,:,1]
        self.array_b = self.array[:,:,2]
        # 2*green/(red+blue)
        self.array_gidx = (2*self.array_g)/(self.array_r+self.array_b+0.01)
        self.h_w_c = self.array.shape

    def crop(self, crp_dim):
        '''
        crop_dim: (left,upper,right,lower)
        '''
        return self.PIL_img.crop(crp_dim)

    def resize(self, resize_dim):
        '''
        resize_dim: (width,height)
        '''
        return self.PIL_img.resize(resize_dim)
    
    def binary_thresh(self, method, value):
        '''
        method: green_channel or green_index
        threshold: the thresh cutoff. 130 for green_channle, 1.12 for green_index
        '''
        if method == 'green_channel':
            _, thresh = cv2.threshold(self.array_g, value, 255, cv2.THRESH_BINARY) # background: white(255), plant: black(0)
            thresh = invert(thresh) # background: black(0), plant: whilte(255)
        elif method == 'green_index':
            _, thresh = cv2.threshold(self.array_gidx, value, 255, cv2.THRESH_BINARY) # background: black(0), plant: whilte(255) 
        thresh = thresh.astype('uint8')
        return thresh # 2d numpy array

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

def toJPG(args):
    '''
    %prog toJPG img1 img2 img3 ...

    convert PNG to JPG using PIL. 
    '''
    p = OptionParser(toJPG.__doc__)
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    for img_fn in args:
        img_out_fn = Path(img_fn).name.replace('.png', '.jpg')
        ProsImage(img_fn).PIL_img.convert('RGB').save(Path(opts.out_dir)/img_out_fn)

def Batch2JPG(args):
    '''
    %prog Batch2JPG in_dir out_dir

    apply toJPG on a large number of images
    '''
    p = OptionParser(Batch2JPG.__doc__)
    p.add_option('--pattern', default='*.png',
                 help="file pattern of png files under the 'dir_in'")
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=Batch2JPG.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    pngs = in_dir_path.glob(opts.pattern)
    cmds = []
    for img_fn in pngs:
        img_fn = str(img_fn).replace(' ', '\ ')
        cmd = "python -m schnablelab.ImageProcessing.base toJPG "\
        f"{img_fn} --out_dir {out_dir}"
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

def Resize(args):
    '''
    %prog Resize img1 img2 img3 ...

    resize image using PIL. 
    If multiple images are provided, same resizing dimension will be applied on all of them
    '''
    p = OptionParser(Resize.__doc__)
    p.add_option('--output_dim', default = '1227,1028',
        help = 'the dimension (width,height) after resizing')
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    p.add_option('--to_jpg', default=False, action='store_true',
        help = 'in save image as jpg format')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    
    dim = ([int(i) for i in opts.output_dim.split(',')])
    for img_fn in args:
        img = ProsImage(img_fn)
        if opts.to_jpg:
            img_out_fn = Path(img.fn).name.replace(f'.{img.format}', '.Rsz.jpg')
            img.resize(dim).convert('RGB').save(Path(opts.out_dir)/img_out_fn)
        else:
            img_out_fn = Path(img.fn).name.replace(f'.{img.format}', f'.Rsz.{img.format}')
            img.resize(dim).save(Path(opts.out_dir)/img_out_fn)

def BatchResize(args):
    '''
    %prog BatchResize in_dir out_dir

    apply BatchResize on a large number of images
    '''
    p = OptionParser(BatchResize.__doc__)
    p.add_option('--pattern', default='*.png',
                 help="file pattern of png files under the 'dir_in'")
    p.add_option('--output_dim', default = '1227,1028',
        help = 'the dimension (width,height) after resizing')
    p.add_option('--to_jpg', default=False, action='store_true',
        help = 'in save image as jpg format')
    p.add_option('--disable_slurm', default=False, action="store_true",
                 help='do not convert commands to slurm jobs')
    p.add_slurm_opts(job_prefix=BatchResize.__name__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    in_dir_path = Path(in_dir)
    pngs = in_dir_path.glob(opts.pattern)
    cmds = []
    for img_fn in pngs:
        img_fn = str(img_fn).replace(' ', '\ ')
        cmd = 'python -m schnablelab.ImageProcessing.base Resize '\
        f'{img_fn} --output_dim {opts.output_dim} --out_dir {out_dir}'
        if opts.to_jpg:
            cmd += ' --to_jpg'
        cmds.append(cmd)
    fn_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    with open(fn_sh, 'w') as f:
        for i in cmds:
            f.write(i+'\n')
    print('check %s for all the commands!'%fn_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

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
        img_fn = str(img_fn).replace(' ', '\ ')
        cmd = "python -m schnablelab.ImageProcessing.base CropObject "\
            f"{img_fn} --out_dir {out_dir} --pad {opts.pad}"
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
        help = 'the dimension (left,upper,right,lower) after cropping. '\
            'using (849,200,1780,1645) for corn images under zoom level 2 ')
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    
    dim = ([int(i) for i in opts.crop_dim.split(',')])
    for img_fn in args:
        img = ProsImage(img_fn)
        img_out_fn = Path(img.fn).name.replace(f'.{img.format}', f'.CrpFrm.{img.format}')
        img.crop(dim).save(Path(opts.out_dir)/img_out_fn)
    
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
        img_fn = str(img_fn).replace(' ', '\ ')
        cmd = "python -m schnablelab.ImageProcessing.base CropFrame "\
            f"{img_fn} --crop_dim {opts.crop_dim} --out_dir {out_dir}"
        cmds.append(cmd)
    cmd_sh = '%s.cmds%s.sh'%(opts.job_prefix, len(cmds))
    pd.DataFrame(cmds).to_csv(cmd_sh, index=False, header=None)
    print('check %s for all the commands!'%cmd_sh)
    if not opts.disable_slurm:
        put2slurm_dict = vars(opts)
        put2slurm(cmds, put2slurm_dict)

def PlantArea(args):
    '''
    %prog PlantArea img1 img2 img3 ...

    estimate plant area 
    '''
    p = OptionParser(PlantArea.__doc__)
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    p.add_option('--frame_size', default = '410,0,2300,1675',
        help = "specify the frame size following left,upper,right,lower. "\
                "For reference: corn image under zoom1: '410,0,2300,1675', "\
                "corn image under zoom2: '849,200,1780,1645'.")
    p.add_option('--thresh_method', default='green_index', choices=('green_index', 'green_channel'),
        help='choose the threshold method')
    p.add_option('--cutoff', type='float', default=1.32,
        help='choose the cutoff of the threshold method. For reference: green_index (1.32), green_channel(130)')
    p.add_option('--out_csv', default='plant_area_summary.csv',
        help='specify csv file including plant pixel area results')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())

    frame_dim = namedtuple('frame_dim', 'left, upper, right, lower')    
    dim = ([int(i) for i in opts.frame_size.split(',')])
    dim = frame_dim._make(dim)

    fns, areas = [], []
    for img_fn in args:
        print(f'Processing {Path(img_fn).name}...')
        fns.append(Path(img_fn).name)
        img = ProsImage(img_fn)
        img_out_fn = Path(img.fn).name.replace(f'.{img.format}', f'.seg.{img.format}')
        arr_thresh = img.binary_thresh(method=opts.thresh_method, value=opts.cutoff)
        arr_thresh[0:dim.upper]=0
        arr_thresh[dim.lower:]=0
        arr_thresh[:, 0:dim.left]=0
        arr_thresh[:, dim.right:]=0
        Image.fromarray(arr_thresh).save(Path(opts.out_dir)/img_out_fn)
        area = (arr_thresh==255).sum()
        areas.append(area)
    df_out = pd.DataFrame(dict(zip(['fn', 'pixel_area'], [fns, areas])))
    df_out['threshold_method'] = opts.thresh_method
    df_out['threshold_cutoff'] = opts.cutoff
    df_out[['fn', 'threshold_method', 'threshold_cutoff', 'pixel_area']].to_csv(Path(opts.out_dir)/opts.out_csv, index=False)
        
if __name__ == "__main__":
    main()