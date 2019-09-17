import numpy as np
from pathlib import Path
from PIL import Image
import cv2
import os
import sys
from sys import argv

def load_3d_image_matrix(hyp_dir, out_fn):
    '''
    Args:
        hyp_dir: string
            filepath of hyperspectral image data
    Returns:
        img_names: 1d array
            File names corresponding to each layer in the img_stack
        img_stack: 3d array
            x,y dims correspond to pixel coordinates for each image
            z dim corresponds to hyperspectral image wavelength.
    '''
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
    if len(sys.argv)==3:
        load_3d_image_matrix(argv[1], argv[2])
    else:
        print('hyp_dir', 'out_fn')
