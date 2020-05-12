# -*- coding: UTF-8 -*-

"""
class and functions to deal with high throughput phenotyping data
"""
import sys
import pandas as pd
from pathlib import Path
from schnablelab.apps.Tools import GenDataFrameFromPath
from schnablelab.apps.base import ActionDispatcher, OptionParser, put2slurm

class ParseProject():
    def __init__(self, prj_dir_name):
        self.prj_dir_name = prj_dir_name
        try:
            df = pd.read_csv('%s.idx.csv'%prj_dir_name)
        except FileNotFoundError:
            print('project index file does not exist, creating one...')
            df = GenDataFrameFromPath(Path(prj_dir_name), pattern='*')
            df['sm'] = df['fn'].apply(lambda x: x.split('_')[1])
            df['date'] = df['fn'].apply(lambda x: x.split('_')[2])
            df['time'] = df['fn'].apply(lambda x: x.split('_')[3])
            df.to_csv('%s.idx.csv'%prj_dir_name, index=False)
        finally:
            self.df = df
        self.sm_counts = self.df['sm'].value_counts().sort_index()
        self.date_counts = self.df['date'].value_counts().sort_index()
        self.SMs = self.df['sm'].unique().values
        self.Dates = self.df['date'].unique().values
    
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
            df = self.df[self.Subsampels]
        elif not samples and dates:
            df = self.df[self.Subdates]
        elif samples and dates:
            df = self.df[self.Subsamples & self.Subdates]
        else:
            df = self.df.copy()
        for fnpath in df['fnpath']:
            yield fnpath.glob('Vis_SV_%s'%angle) if angle else fnpath.glob('Vis_*')
            
def main():
    actions = (
        ('ExtractRGBs', 'extract images from project folder'),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

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
    opts.samples = opts.samples.split(',') if opts.samples else opts.samples
    opts.dates = opts.dates.split(',') if opts.dates else opts.dates

    prj = ParseProject(project_folder)
    RGBs = prj.RGB(samples=opts.samples, dates=opts.dates, angle=opts.angle)
    for rgb in RGBs:
        print(rgb)

if __name__ == '__main__':
    main()