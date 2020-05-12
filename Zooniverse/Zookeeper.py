'''
schnablelab CLI tool
Calls Zookeeper class
'''
import sys
import math
from pandas as pd
from shutil import copyfile
from schnablelab.apps.Tools import GenDataFrameFromPath
from schnablelab.apps.base import ActionDispatcher, OptionParser, cutlist

def main():
    actions = (
        ('divide', 'divide a large number of images to sevearl smaller sets'),
        ('upload', 'load images to zooniverse'),
        ('export', 'Get annotation and other exports'),
        ('manifest', 'Generate a manifest for zooniverse subject set upload')
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

def divide(args):
    '''
    %prog divide input_dir output_dir_prefix
    '''
    p = OptionParser(divide.__doc__)
    p.add_option('--nimgs_per_folder', type='int', default=900,
                 help='~ number of images in each smaller folder.')
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    input_dir, out_prefix, args

    df = GenDataFrameFromPath(Path(input_dir), pattern='*.png')
    n_folders = math.ceil(df.shape[0]/opts.nimgs_per_folder)
    n = 0
    for _, grp in cutlist(df['fnpath'].values, n_folders):
        n += 1
        output_folder = Path('%s_'%(out_prefix,n))
        print(output_folder)
        if not output_folder.exists():
            output_folder.mkdir()
        for i in grp:
            copyfile(i, out_dir/i.name)

def upload(args):
    '''
    %prog upload imgdir projid

    - imgdir: Path to directory of the images to be uploaded
    - projid: Zooniverse project id (4 - 5 digit number)

    DESC:
        Uploads images from the image directory to zooniverse
        project. If there is no manifest will generate one.
    '''

    from schnablelab.Zooniverse.Zootils import upload as load

    p = OptionParser(upload.__doc__)
    p.add_option('-s', '--subject', default=False,
                 help='Designate a subject set id.')
    p.add_option('-q', '--quiet', action='store_true', default=False,
                 help='Silences output when uploading images to zooniverse.')
    p.add_option('-x', '--extension', default=False,
                 help='Specify the extension of the image files to be uploaded.')
    '''
    p.add_option('-c', '--convert', action='store_true', default=False,
                 help="Compress and convert files to jpg for faster load times"
                 + " on zooniverse.\n"
                 + " Command: magick -strip -interlace Plane -quality 85%"
                 + " -format jpg <img_directory>/<filename>.png")
    '''

    opts, args = p.parse_args(args)

    if len(args) != 2:
        p.print_help()
        exit(False)

    imgdir, projid = args

    load(imgdir, projid, opts)

    return True


def export(args):
    '''
    %prog export proj_id outfile

    - proj_id: The project id of the zooniverse project

    DESC: Fetches an export from the specified zooniverse project id.
    '''

    from schnablelab.Zooniverse.Zootils import export as exp

    p = OptionParser(export.__doc__)
    p.add_option('-t', '--type', default='classifications',
                 help='Specify the type of export')

    opts, args = p.parse_args(args)

    if len(args) != 2:
        exit(not p.print_help())

    projid, outfile = args

    exp(projid, outfile, opts)

    return True


def manifest(args):
    '''
    %prog manifest image_dir

    - img_dir: The image directory in which to generate the manifest.

    DESC: Generates a manifest inside the specified image directory.
    '''
    from schnablelab.Zooniverse.Zootils import manifest as mani

    p = OptionParser(manifest.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        exit(not p.print_help())

    imgdir = args[0]

    mani(imgdir)

    return True


if __name__ == '__main__':
    main()
