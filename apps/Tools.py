from __future__ import print_function
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import linregress
from collections import defaultdict
from sklearn.metrics import mean_squared_error

def print_json(data):
    print(json.dumps(data, indent=2))

def bin_markers(df, diff=0, missing_value='-'):
    """
    merge consecutive markers with same genotypes
    return slelected row index
    Examples:
    """
    df = df.replace(missing_value, np.nan)
    first_row = df.iloc[0,:]
    temp = [df.index[0]] # save temp index
    pre_row = first_row 
    df_rest = df.iloc[1:,:]
    result_ids = []
    for idx, row in df_rest.iterrows():
        df_tmp = pd.concat([pre_row, row], axis=1).dropna()
        diff_num = (df_tmp.iloc[:,0] != df_tmp.iloc[:,1]).sum()
        if diff_num <= diff:
            temp.append(idx)
        else:
            if len(temp) > 1:
                result_ids.append(temp)
            else:
                result_ids.append([idx])
            temp = [idx]
        pre_row = row
    if result_ids[0][0] != df.index[0]:
        result_ids.insert(0, [df.index[0]])

    results = []
    represent_idx, block_idx = [], []
    for index in result_ids:
        if len(index) > 1:
            df_tmp = df.loc[index, :]
            good_idx = df_tmp.isnull().sum(axis=1).idxmin()
            results.append(good_idx)
            represent_idx.append(good_idx)
            block_idx.append(index)
        else:
            results.append(index[0])
    return represent_idx, block_idx, results




def sort_merge_sort(arrays):
    """
    get redundant lists by merging lists with overlaping region.
    Example:
    >>> a = [[1,3], [3, 5], [6,10], [7, 9], [11,15], [11,12],[16,30]]
    >>> sort_merge_sort(a)
    >>> [array([1, 3, 5]), array([ 6,  7,  9, 10]), array([11, 12, 15]), [16, 30]]
    """
    val_start = [i[0] for i in arrays]
    val_end = [i[-1] for i in arrays]
    df = pd.DataFrame(dict(zip(['array', 'val_start', 'val_end'], [arrays, val_start, val_end]))).sort_values(['val_start', 'val_end']).reset_index(drop=True)
    first_arr = df.loc[0, 'array']
    temp = first_arr
    pre_arr = first_arr
    results = []
    for arr in df.loc[1:, 'array']:
        if arr[0] <= pre_arr[-1]:
            temp.extend(arr)
        else:
            if len(temp) == len(pre_arr):
                results.append(pre_arr)
            else:
                temp_sorted_unique = pd.Series(temp).sort_values().unique()
                results.append(temp_sorted_unique)
            temp = arr
        pre_arr = arr
    if len(temp) == len(pre_arr):
        results.append(pre_arr)
    else:
        temp_sorted_unique = pd.Series(temp).sort_values().unique()
        results.append(temp_sorted_unique)
    return results

def get_blocks(np_1d_array, dist=150, block_size=2):
    """
    group values to a block with specified distance
    Examples:
    >>> a = np.array([1,2,4,10,12,13,15])
    >>> test(a, dist=1)
    [[1, 2], [12, 13]]
    >>> test(a, dist=2)
    [[1, 2, 4], [10, 12, 13, 15]]
    """
    first_val = np_1d_array[0]
    temp = [first_val] # save temp blocks
    pre_val = first_val 
    results = []
    for val in np_1d_array[1:]:
        if (val - pre_val) <= dist:
            temp.append(val)
        else:
            if len(temp) >= block_size:
                results.append(temp)
            temp = [val]
        pre_val = val
    if len(temp) >= block_size:
        results.append(temp)
    return results

def random_alternative(lens, values=[0,2]):
    """
    return a numpy array with alternating interger values
    """
    v1, v2 = values
    st_value = np.random.choice(values)
    alternative_value = v1 if st_value == v2 else v2
    a = np.empty((lens,))
    a[::2] = st_value
    a[1::2] = alternative_value
    return a.astype('int')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def GenDataFrameFromPath(path, pattern='*.png'):
    """
    generate a dataframe for all file in a dir with the specific pattern of file name.
    use: GenDataFrameFromPath(path, pattern='*.png')
    """
    fnpaths = list(path.glob(pattern))
    df = pd.DataFrame(dict(zip(['fnpath'], [fnpaths])))
    df['dir'] = df['fnpath'].apply(lambda x: x.parent)
    df['fn'] = df['fnpath'].apply(lambda x: x.name)
    return df

def ConciseVcf(fn):
    """
    concise the vcf file by remove the header, useless columns and simplfied genotype
    ConciseVcf(fn)
    """
    n = 0
    f = open(fn)
    for i in f:
        if i.startswith('##'):
            n += 1
        else:
            break
    df = pd.read_csv(fn, header=n, delim_whitespace=True)
    df = df.drop(['INFO', 'FORMAT', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER'], axis=1)
    for idx in df.columns[2:]:
        df[idx] = df[idx].map(lambda x: x.split(':')[0])
    df = df.replace(['0/0', '0/1', '1/0', '1/1', './.'], [0, 1, 1, 2, 9])
    return df

def getChunk(fn, ignore=1):
    '''ignore: rows starts with pound sign'''
    df0_chr = defaultdict(int)
    chr_order = []
    with open(fn) as f:
        for dash_line in range(ignore):
            f.readline()
        for i in f:
            j = i.split()[0].split('-')[0]
            df0_chr[j] += 1
            if j in chr_order:
                pass
            else:
                chr_order.append(j)
    if len(chr_order) != len(set(chr_order)):
        sys.exit('Please check your marker name and sort them by chr name.')
    return chr_order, df0_chr
            
class SimpleStats(object):
    """
    This class will do the simple statistics on two series objecjts.
    a) linear regressoin: slope, intercept, r^2, p_value
    b) mean, std of the difference and absolute differnece
    c) MSE (mean squared error) and RMSE (root mean squared error)
    d) agreement
    e) plot the regreesion figure and the difference distribution figure
    """
    def __init__(self, series1, series2):
        self.s1 = series1
        self.s2 = series2
        self.length = series1.shape[0]
        self.diff = series1 - series2
        self.absdiff = (series1 - series2).abs()

    def regression(self):
        slope, intercept, r_value, p_value, __ = linregress(self.s1, self.s2)
        return slope, intercept, r_value**2, p_value

    def mean_std_diff(self):
        mean, std = self.diff.mean(), self.diff.std()
        return mean, std

    def mean_std_absdiff(self):
        abs_mean, abs_std = self.absdiff.mean(), self.absdiff.std()
        return abs_mean, abs_std

    def mse(self):
        mse = mean_squared_error(self.s1, self.s2)
        return mse

    def rmse(self):
        rmse = mean_squared_error(self.s1, self.s2)**0.5
        return rmse
    
    def agreement(self, cutoff):
        return (self.absdiff<=float(cutoff)).sum()/self.length

































