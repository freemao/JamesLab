# 7/16/18
# chenyong 
# prediction

"""
make predictions using trained model
"""
import os
import math
import os.path as op
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from PIL import Image
from schnablelab.apps.base import cutlist, ActionDispatcher, OptionParser, glob
from schnablelab.apps.headers import Slurm_header, Slurm_gpu_header
from schnablelab.apps.natsort import natsorted
from glob import glob
from pathlib import Path
from subprocess import run

def main():
    actions = (
        ('Plot', 'plot training model history'),
        ('Predict', 'using trained neural network to make prediction'),
        ('Predict_Rmodels', 'using trained model to make prediction'),
        ('Predict_Rmodels_slurms', 'generate slurm  jobs for predict_Rmodels'),
        ('PredictSlurmCPU', 'generate slurm CPU job of prediction'),
        ('PredictSlurmGPU', 'generate slurm GPU job of prediction'),
        ('Imgs2Arrs', 'convert hyperspectral images under a dir to a numpy array object'),
        ('Imgs2ArrsBatch', 'generate slurm job convert HyperImages to NPY for all hyperspectral image dirs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def Predict_Rmodels(args):
    """
    %prog Predict_Rmodels model_name input_2d_npy/csv corresponding_3d_npy
    using R models to make prediciton on the whole hyper image
    """
    p = OptionParser(Predict_Rmodels.__doc__)
    opts, args = p.parse_args(args)
    if len(args) != 3:
        sys.exit(not p.print_help())
    model_name, npy2d, npy3d, = args
    h,_,_ = np.load(npy3d).shape
    Rfile = op.abspath(op.dirname(__file__)) + '/R_models/%s_Predict.R'%model_name
    cmd0 = 'ml R'
    cmd1 = 'R CMD BATCH -%s -%s %s'%(npy2d, h, Rfile)
    print(cmd0)
    print(cmd1)
    run(cmd1, shell=True)

def Predict_Rmodels_slurms(args):
    '''
    %prog Predict_Rmodels_slurms input_2dcsv/npy_dir corresponding_3dnpy_dir
    
    generate slurm file for Predict_Rmodels
    '''
    p = OptionParser(Predict_Rmodels_slurms.__doc__)
    p.add_option("--model", default='LDA', choices=('LDA', 'MLR', 'SVM', 'PLSDA', 'RF', 'QDA', 'LASSO'),
        help="input 2d array format")
    p.add_option("--input_format", default='npy', choices=('csv', 'npy'),
        help="input 2d array format")
    p.add_option('--pattern', default='*.2d.npy',choices=('*.2d.csv', '*.2d.npy'),
        help='2d npy/csv file patterns in dir')
    p.add_option("--n", default=10,type='int',
        help="number of commands in each slurm file")
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())
    in_dir, npy3d_dir, = args

    dir_path = Path(in_dir)
    csvs = list(dir_path.glob(opts.pattern))
    csvs_n = len(csvs)

    npy_dir = Path(npy_dir)
    print('%s 2d array files found...'%csvs_n)
    for i in range(math.ceil(csvs_n/opts.n)):
        batch_csvs = csvs[i*opts.n: (i+1)*opts.n]
        print('batch%s'%i, len(batch_csvs))
        cmd = ''
        for csv in batch_csvs:
            npy3d_fn = csv.name.replace('.2d.csv', '.npy') if opts.input_format=='csv' else csv.name.replace('.2d.npy', '.npy')
            npy3d = npy3d_dir/npy3d_fn
            cmd += 'python -m schnablelab.CNN.Predict_snn Predict_Rmodels %s %s %s\n'%(opts.model, csv, npy3d)
        prefix = '%s_Predict_batch%s'%(opts.model, i)
        header = Slurm_header%(10, 10000, prefix, prefix, prefix)
        header += 'conda activate MCY\n'
        header += 'module load R\n'
        header += cmd
        with open('%s.slurm'%prefix, 'w') as f:
            f.write(header)

def Imgs2Arrs(args):
    '''
    %prog hyp_dir(filepath of hyperspectral image data) 
    Returns: numpy array object with shape [x*y, z].
        x,y dims correspond to pixel coordinates for each image
        z dim corresponds to hyperspectral image wavelength.
    '''
    import cv2
    
    p = OptionParser(Imgs2Arrs.__doc__)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    mydir, = args
    imgs = [i for i in os.listdir(mydir) if i.endswith('png')]
    sorted_imgs = sorted(imgs, key=lambda x: int(x.split('_')[0]))
    all_arrs = []
    for i in sorted_imgs[2:]:
        print(i)
        #img = cv2.imread('%s/%s'%(mydir, i), cv2.IMREAD_GRAYSCALE)
        img = np.array(Image.open('%s/%s'%(mydir, i)).convert('L'))
        print(img.shape)
        all_arrs.append(img)
    arrs = np.stack(all_arrs, axis=2)
    np.save('%s.npy'%mydir, arrs)

def Imgs2ArrsBatch(args):
    """
    %prog HyperDirPattern("CM*")
    generate img2arr jobs for all hyperspectral image dirs
    """
    p = OptionParser(Imgs2ArrsBatch.__doc__)
    p.set_slurm_opts()
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    pattern, = args
    all_dirs = [i for i in glob(pattern) if os.path.isdir(i)]
    for i in all_dirs:
        cmd = 'python -m schnablelab.CNN.Predict Imgs2Arrs %s\n'%i
        jobname = i+'.img2npy'
        header = Slurm_header%(opts.time, opts.memory, jobname, jobname, jobname)
        #header += "ml anaconda\nsource activate MCY\n"
        header += cmd
        jobfile = open('%s.img2arr.slurm'%i, 'w')
        jobfile.write(header)
        jobfile.close()
        print('slurm job for %s has been generated.'%i)

def Predict(args):
    """
    %prog model_name npy_pattern('CM*.npy')
    using your trained model to make predictions on selected npy (2d or 3d) files.
    The pred_data is a numpy array object which has the same number of columns as the training data.
    """
    from keras.models import load_model
    import scipy.misc as sm
    p = OptionParser(Predict.__doc__)
    p.add_option('--range', default='all',
        help = "specify the range of the testing images, hcc job range style")
    p.add_option('--opf', default='infer',
        help = "specify the prefix of the output file names")
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())
    model, npy_pattern = args
    opf = model.split('/')[-1].split('.')[0] if opts.opf == 'infer' else opts.opf

    npys = glob(npy_pattern)
    if opts.range != 'all':
        start = int(opts.range.split('-')[0])
        end = int(opts.range.split('-')[1])
        npys = npys[start:end]
    print('%s npys will be predicted this time.'%len(npys))

    my_model = load_model(model)
    for npy in npys:
        print(npy)
        test_npy = np.load(npy)
        npy_shape = test_npy.shape
        np_dim = len(npy_shape)
        test_npy_2d = test_npy.reshape(npy_shape[0]*npy_shape[1], npy_shape[2]) if np_dim==3 else test_npy
        print('testing data shape:', test_npy_2d.shape)
        pre_prob = my_model.predict(test_npy_2d/255)
        predictions = pre_prob.argmax(axis=1) # this is a numpy array

        if np_dim == 3:
            predictions = predictions.reshape(npy_shape[0], npy_shape[1])
            df = pd.DataFrame(predictions)
            df1 = df.replace(0, 255).replace(1, 127).replace(2, 253).replace(3, 190)#0: background; 1: leaf; 2: stem; 3: panicle
            df2 = df.replace(0, 255).replace(1, 201).replace(2, 192).replace(3, 174)
            df3 = df.replace(0, 255).replace(1, 127).replace(2, 134).replace(3, 212) 
            arr = np.stack([df1.values, df2.values, df3.values], axis=2)
            opt = npy.split('/')[-1].split('.npy')[0]+'.prd'
            sm.imsave('%s.%s.png'%(opf,opt), arr)
        elif np_dim == 2:
            opt = npy.split('/')[-1].split('.npy')[0]+'.prd'
            np.savetxt('%s.%s.csv'%(opf, opt), predictions)
        else:
            sys.exit('either 2 or 3 dim numpy array!')
        print('Done!')

def PredictSlurmCPU(args):
    """
    %prog model_name npyPattern("CM*.npy") job_n
    generate prediction CPU jobs for all npy files
    """
    p = OptionParser(PredictSlurmCPU.__doc__)
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    #print(args)
    mn, npy_pattern, jobn, = args
    if opts.prefix == 'myjob':
        print('specify job name prefix!') 
        sys.exit()

    npys = glob(npy_pattern)
    print(len(npys))
    grps = cutlist(npys, int(jobn))
    for gn, grp in grps:
        st, ed = gn.split('-')
        ed = int(ed)+1
        gn = '%s-%s'%(st, ed)
        cmd = "python -m schnablelab.CNN.Predict_snn Predict %s '%s'\n"%(mn, npy_pattern)
        opt = '%s.%s'%(opts.prefix, gn)
        header = Slurm_header%(opts.time, opt, opt, opt, opt)
        header += "ml anaconda\nsource activate Py3KerasTensorCPU\n"
        header += cmd
        with open('%s.cpu.slurm'%opt, 'w') as f:
            f.write(header)
        print('%s.cpu.slurm prediction CPU job file generated!'%opt)

def PredictSlurmGPU(args):
    """
    %prog model_name npyPattern("CM*.npy") job_n
    generate prediction GPU jobs for all npy files
    """
    p = OptionParser(PredictSlurmGPU.__doc__)
    p.set_slurm_opts(jn=True)
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    mn, npy_pattern, jobn, = args
    if opts.prefix == 'myjob':
        print('specify job name prefix!') 
        sys.exit()

    npys = glob(npy_pattern)
    print(len(npys))
    grps = cutlist(npys, int(jobn))
    for gn, grp in grps:
        st, ed = gn.split('-')
        ed = int(ed)+1
        gn = '%s-%s'%(st, ed)
        cmd = "python -m schnablelab.CNN.Predict_snn Predict %s '%s'\n"%(mn, npy_pattern)
        opt = '%s.%s'%(opts.prefix, gn)
        header = Slurm_gpu_header%(opts.time, opts.memory, opt, opt, opt)
        header += "ml anaconda\nsource activate MCY\n"
        header += cmd
        with open('%s.gpu.slurm'%opt, 'w') as f:
            f.write(header)
        print('%s.gpu.slurm prediction GPU job file generated!'%opt)

def Plot(args): 
    """
    %prog dir
    plot training process
    You can load the dict back using pickle.load(open('*.p', 'rb'))
    """

    p = OptionParser(Plot.__doc__)
    p.add_option("--pattern", default="History_*.p",
        help="specify the pattern of your pickle object file, remember to add quotes [default: %default]")
    p.set_slurm_opts()
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    mydir, = args
    pickles = glob('%s/%s'%(mydir, opts.pattern)) 
    print('total %s pickle objects.'%len(pickles))
    #print(pickles)
    for p in pickles:
        fs, es = opts.pattern.split('*')
        fn = p.split(fs)[-1].split(es)[0]
        myp = pickle.load(open(p, 'rb'))
        
        mpl.rcParams['figure.figsize']=[7.5, 3.25]
        fig, axes = plt.subplots(nrows=1, ncols=2)
        
        # summarize history for accuracy
        ax1 = axes[0]
        ax1.plot(myp['acc'])
        ax1.plot(myp['val_acc'])
        ax1.set_title('model accuracy')
        ax1.set_ylabel('accuracy')
        ax1.set_xlabel('epoch')
        ax1.set_ylim(0,1.01)
        ax1.legend(['train', 'validation'], loc='lower right')
        max_acc = max(myp['val_acc'])
	# summarize history for loss
        ax2 = axes[1]
        ax2.plot(myp['loss'])
        ax2.plot(myp['val_loss'])
        ax2.set_title('model loss')
        ax2.set_ylabel('loss')
        ax2.set_xlabel('epoch')
        ax2.legend(['train', 'validation'], loc='upper right')
        plt.tight_layout()
        plt.savefig('%s_%s.png'%(max_acc,fn))    
        plt.clf()

if __name__ == "__main__":
    main()
