# -*- coding: UTF-8 -*-

"""
class and functions to deal with high throughput phenotyping data
"""
import sys
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from shutil import copyfile
from schnablelab.apps.Tools import GenDataFrameFromPath
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

class ParseProject():
    def __init__(self, prj_dir_name):
        self.prj_dir_name = prj_dir_name
        try:
            df = pd.read_csv('%s.idx.csv'%prj_dir_name)
            df['fnpath'] = df['fnpath'].apply(lambda x: Path(x))
        except FileNotFoundError:
            print('project index file does not exist, creating one...')
            df = GenDataFrameFromPath(Path(prj_dir_name), pattern='*')
            df['sm'] = df['fn'].apply(lambda x: x.split('_')[1])
            df['date'] = df['fn'].apply(lambda x: x.split('_')[2])
            df['time'] = df['fn'].apply(lambda x: x.split('_')[3])
            df = df.sort_values(['sm', 'date', 'time']).reset_index(drop=True)
            df.to_csv('%s.idx.csv'%prj_dir_name, index=False)
        finally:
            self.df = df
        self.sm_counts = self.df['sm'].value_counts().sort_index()
        self.date_counts = self.df['date'].value_counts().sort_index()
        self.SMs = self.df['sm'].unique()
        self.Dates = self.df['date'].unique()
    
    def Subsamples(self, samples):
        '''
        samples (list): list of samples
        '''
        for sm in samples:
            if not sm in self.SMs:
                sys.exit('%s not in the sample list'%sm)
        cond = self.df['sm'].isin(samples)
        return cond
    
    def Subdates(self, dates):
        '''
        dates (list): list of dates
        '''
        for date in dates:
            if not date in self.Dates:
                sys.exit('%s not in the date list'%date)
        cond = self.df['date'].isin(dates)
        return cond

    def RGB(self, samples=None, dates=None, angle=None):
        if samples and not dates:
            df = self.df[self.Subsamples(samples)]
        elif not samples and dates:
            df = self.df[self.Subdates(dates)]
        elif samples and dates:
            df = self.df[self.Subsamples(samples) & self.Subdates(dates)]
        else:
            df = self.df.copy()
        pbar = tqdm(df.iterrows(), total=df.shape[0])
        for _,row in pbar:
            sm, d, hms = row['sm'], row['date'], row['time']
            results = row['fnpath'].glob('Vis_SV_%s'%angle) if angle else row['fnpath'].glob('Vis_*')
            pbar.set_description('extracting %s %s %s...'%(sm, d, hms))
            yield sm, d, hms, results
            
def main():
    actions = (
        ('ExtractRGBs', 'extract images from project folder'),
        ('Info', 'summary of image data under the project folder'),
        ('List', 'list specified image folders'),

    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def List(args):
    '''
    %prog List project_folder
    
    list specified image folders
    '''
    p = OptionParser(List.__doc__)
    p.add_option('--samples',
        help = 'specify samples (comma separated without space)')
    p.add_option('--dates',
        help = 'specify dates (comma separated without space)')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    project_folder, = args

    prj = ParseProject(project_folder)
    if opts.samples and not opts.dates:
        samples = opts.samples.split(',')
        cond = prj.Subsamples(samples)
    elif not opts.samples and opts.dates:
        dates = opts.dates.split(',')
        cond = prj.Subdates(dates)
    elif not opts.samples and not opts.dates:
        print('Specify either samples or dates for showing!')
    else:
        print('provide either samples or dates for showing')
    print(prj.df[cond][['fn', 'sm', 'date', 'time']])

def Info(args):
    '''
    %prog Info project_folder

    Show summary of images under project_folder
    '''
    p = OptionParser(Info.__doc__)
    opts, args = p.parse_args(args)
    
    if len(args) == 0:
        sys.exit(not p.print_help())
    project_folder, = args

    prj = ParseProject(project_folder)
    print('Summary of samples:')
    print(prj.sm_counts, '\n')
    print('Summary of dates:')
    print(prj.date_counts, '\n')
    print('Angles for RGB images:')
    for angle in prj.df.loc[0,'fnpath'].glob('Vis_*'):
        print(angle.name)

def ExtractRGBs(args):
    '''
    %prog ExtractRGBs project_folder

    extract RGB images from project folder
    '''
    p = OptionParser(ExtractRGBs.__doc__)
    p.add_option('--out_dir', default='.',
        help = 'specify the output image directory')
    p.add_option('--samples',
        help = 'samples')
    p.add_option('--dates',
        help = 'dates')
    p.add_option('--angle',
        help = 'RBG viewing angle')
    opts, args = p.parse_args(args)
    
    if len(args) == 0:
        sys.exit(not p.print_help())
    project_folder, = args

    out_dir = Path(opts.out_dir)
    if not out_dir.exists():
        print('%s does not exist, creating..'%out_dir)
        out_dir.mkdir()

    opts.samples = opts.samples.split(',') if opts.samples else opts.samples
    opts.dates = opts.dates.split(',') if opts.dates else opts.dates

    prj = ParseProject(project_folder)
    for sm, d, hms, RGBs in prj.RGB(samples=opts.samples, dates=opts.dates, angle=opts.angle):
        for rgb in RGBs:
            out_fn = '%s_%s_%s_%s.png'%(sm, d, hms, rgb.name)
            copyfile(rgb/'0_0_0.png', out_dir/out_fn)
        
if __name__ == '__main__':
    main()
