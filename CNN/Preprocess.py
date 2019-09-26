import numpy as np
import math
from pathlib import Path
from PIL import Image
import cv2
import os
import sys
from sys import argv
from schnablelab.apps.natsort import natsorted
from schnablelab.apps.headers import Slurm_header
from schnablelab.apps.base import ActionDispatcher, OptionParser, glob, iglob

def main():
    actions = (
        ('crop_png', 'crop png images'),
        ('crop_png_slurms', 'gen slurm jobs for crop_png'),
        ('hyp2arr', 'convert hyperspectral images to a numpy array'),
        ('hyp2arr_slurms', 'gen slurm jobs for hyp2arr'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def crop_png(args):
    '''
    %prog crop_png fn_in fn_out coordinates(left,upper,right,lower)
    '''
    p = OptionParser(crop_png.__doc__)
    opts, args = p.parse_args(args)
    fn_in, fn_out, crds, = args
    crp_tuple = ([int(i) for i in crds.split(',')])
    img = Image.open(fn_in).crop(crp_tuple)
    img.save(fn_out, 'PNG')

def crop_png_slurms(args):
    '''
    %prog crop_png_slurms in_dir out_dir
    
    generate slurm file for crop_png
    '''
    p = OptionParser(crop_png_slurms.__doc__)
    p.add_option('--pattern', default='*.png',
        help='image file patterns in_dir')
    p.add_option("--crop_crds", default='760,0,1891,1671',
        help="left,upper,right,lower pixel coordinate, separated by comma")
    p.add_option("--n", default=10,type='int',
        help="number of commands in each slurm file")

    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        print('%s does not exist, create...')
        out_path.mkdir()
    dir_path = Path(in_dir)
    pngs = list(dir_path.glob(opts.pattern))
    pngs_n = len(pngs)
    print('%s images found...'%pngs_n)
    
    for i in range(math.ceil(pngs_n/opts.n)):
        batch_pngs = pngs[i*opts.n: (i+1)*opts.n]
        print('batch%s'%i, len(batch_pngs))
        cmd = ''
        for png in batch_pngs:
            out_png = png.name.replace('.png', '.crp.png')
            out_png_path = out_path/out_png
            cmd += 'python -m schnablelab.CNN.Preprocess crop_png %s %s %s\n'%(png, out_png_path, opts.crop_crds)
        prefix = 'crp_batch%s'%i
        header = Slurm_header%(10, 1000, prefix, prefix, prefix)
        header += 'conda activate MCY\n'
        header += cmd
        with open('%s_crop_%s.crop.slurm'%(in_dir, i), 'w') as f:
            f.write(header)

def hyp2arr_slurms(args):
    '''
    %prog hyp2arr_slurms in_dir out_dir
    
    generate hyp2arr slurm jobs for all folders under specified dir
    '''
    p = OptionParser(hyp2arr_slurms.__doc__)
    p.add_option('--pattern', default='*',
                 help='hyper dir pattern for folders under dir')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    in_dir, out_dir, = args
    out_path = Path(out_dir)
    if not out_path.exists():
        sys.exit('%s does not exist...')
    dir_path = Path(in_dir)
    folders = list(dir_path.glob(opts.pattern))
    num_arrs = len(folders)
    print('%s hyper folders found'%num_arrs)
    for hyp_dir in folders:
        in_dir = str(hyp_dir/'Hyp_SV_90')
        out_fn = hyp_dir.name.replace(' ', '_')
        out_fn_path = out_path/out_fn
        cmd = 'python -m schnablelab.CNN.Preprocess hyp2arr %s %s'%(in_dir, out_fn_path)
        print(cmd)
        header = Slurm_header%(10, 5000, out_fn, out_fn, out_fn)
        header += 'conda activate MCY\n'
        header += cmd
        with open('%s.hyp2arr.slurm'%out_fn, 'w') as f:
            f.write(header)

def hyp2arr(args):
    '''
    %prog hyp2arr hyp_dir out_fn

    convert hyperspectral images to numpy array
    '''
    p = OptionParser(hyp2arr.__doc__)
    opts, args = p.parse_args(args)
    if len(args)==0:
        sys.exit(not p.print_help())
    hyp_dir, out_fn, = args

    discard_imgs = ['0_0_0.png', '1_0_0.png']
    dir_path = Path(hyp_dir)
    if not dir_path.exists():
        sys.exit('%s does not exist!'%hyp_dir)
    imgs = list(dir_path.glob('*.png'))
    imgs = sorted(imgs, key=lambda x: int(x.name.split('_')[0]))
    num_imgs = len(imgs)
    print('%s images found.'%num_imgs)
    img_arrs = []
    for i in imgs:
        if not i.name in discard_imgs:
            arr = cv2.imread(str(i), cv2.IMREAD_GRAYSCALE)
            img_arrs.append(arr)
    img_array = np.stack(img_arrs, axis=2)
    print(img_array.shape)
    np.save(out_fn, img_array)

if __name__=='__main__':
    main()
