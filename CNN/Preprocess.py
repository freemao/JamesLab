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
        ('three2two', 'convert 3d npy to 2d'),
        ('three2two_slurms', 'gen slurm jobs for three2two'),
        ('crop_png', 'crop png images'),
        ('crop_png_slurms', 'gen slurm jobs for crop_png'),
        ('resize_png', 'resize png images'),
        ('resize_png_slurms', 'gen slurm jobs for resize_png'),
        ('hyp2arr', 'convert hyperspectral images to a numpy array'),
        ('hyp2arr_slurms', 'gen slurm jobs for hyp2arr'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def three2two_slurms(args):
    '''
    %prog three2two_slurms in_dir out_dir
    
    generate slurm file for three2two
    '''
    p = OptionParser(three2two_slurms.__doc__)
    p.add_option('--pattern', default='*.npy',
        help='npy file patterns in_dir')
    p.add_option('--crops',
        help='the coordinates for croping, follow left,upper,right,lower format. 1,80,320,479')
    p.add_option("--out_format", default='npy', choices=('npy', 'csv'),
        help="choose the output format")
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
    npys = list(dir_path.glob(opts.pattern))
    npys_n = len(npys)
    print('%s npy files found...'%npys_n)
    for i in range(math.ceil(npys_n/opts.n)):
        batch_npys = npys[i*opts.n: (i+1)*opts.n]
        print('batch%s'%i, len(batch_npys))
        cmd = ''
        for npy in batch_npys:
            out_prefix = npy.name.replace('.npy', '')
            out_npy_path = out_path/out_prefix
            if opts.crops:
                cmd += 'python -m schnablelab.CNN.Preprocess three2two %s %s --crops %s --format %s\n'%(npy, out_npy_path, opts.crops, opts.out_format)
            else:
                cmd += 'python -m schnablelab.CNN.Preprocess three2two %s %s --format %s\n'%(npy, out_npy_path, opts.out_format)
        prefix = 'three2two_batch%s'%i
        header = Slurm_header%(10, 10000, prefix, prefix, prefix)
        header += 'conda activate MCY\n'
        header += cmd
        with open('three2two_%s.slurm'%i, 'w') as f:
            f.write(header)

def three2two(args):
    '''
    %prog three2two fn_in out_prefix 

    convert 3d npy to 2d
    '''
    p = OptionParser(three2two.__doc__)
    p.add_option('--crops',
        help='the coordinates for croping, follow left,upper,right,lower format. 1,80,320,479')
    p.add_option("--format", default='npy', choices=('npy', 'csv'),
        help="choose the output format")
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    fn_in, out_prefix, = args
    npy = np.load(fn_in)
    if opts.crops:
        left, up, right, down = opts.crops.split(',')
        npy = npy[int(up):int(down),int(left):int(right),:]
    h,w,d = npy.shape
    print(h, w, d)
    npy_2d = npy.reshape(h*w, d)
    if opts.format=='csv':
        out_fn = "%s.2d.csv"%out_prefix
        np.savetxt(out_fn, npy_2d, delimiter=",")
    else:
        out_fn = "%s.2d.npy"%out_prefix
        np.save(out_fn, npy_2d.astype(np.float64))
    print('Done!')

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

def resize_png(args):
    '''
    %prog resize_png fn_in fn_out dimension(withd,height)
    '''
    p = OptionParser(resize_png.__doc__)
    opts, args = p.parse_args(args)
    fn_in, fn_out, crds, = args
    crp_tuple = ([int(i) for i in crds.split(',')])
    img = Image.open(fn_in).resize(crp_tuple)
    img.save(fn_out, 'PNG')

def resize_png_slurms(args):
    '''
    %prog resize_png_slurms in_dir out_dir
    
    generate slurm file for resize_png
    '''
    p = OptionParser(resize_png_slurms.__doc__)
    p.add_option('--pattern', default='*.png',
        help='image file patterns in_dir')
    p.add_option("--resize_dim", default='250,250',
        help="target width and height separated by comma")
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
            out_png = png.name.replace('.png', '.resize.png')
            out_png_path = out_path/out_png
            cmd += 'python -m schnablelab.CNN.Preprocess resize_png %s %s %s\n'%(png, out_png_path, opts.resize_dim)
        prefix = 'resize_batch%s'%i
        header = Slurm_header%(10, 1000, prefix, prefix, prefix)
        header += 'conda activate MCY\n'
        header += cmd
        with open('%s_resize_%s.crop.slurm'%(in_dir.split('/')[-1], i), 'w') as f:
            f.write(header)

def crop_png_slurms(args):
    '''
    %prog crop_png_slurms in_dir out_dir
    
    generate slurm file for crop_png
    '''
    p = OptionParser(crop_png_slurms.__doc__)
    p.add_option('--pattern', default='*.png',
        help='image file patterns in_dir')
    p.add_option("--crop_crds", default='1000,2800,3500,5455',
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
        with open('%s_crop_%s.crop.slurm'%(in_dir.split('/')[-1], i), 'w') as f:
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
            print(i)
            arr = cv2.imread(str(i), cv2.IMREAD_GRAYSCALE)
            print(i.name, arr.shape)
            img_arrs.append(arr)
    img_array = np.stack(img_arrs, axis=2)
    print(img_array.shape)
    np.save(out_fn, img_array)

if __name__=='__main__':
    main()
