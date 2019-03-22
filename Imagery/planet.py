#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
search, stat, and download data from Planet through planet api endpoints
"""
import os
import sys
import json
import time
import logging
import requests
import geojsonio
import os.path as op
import pandas as pd
import numpy as np
from pathlib import Path
from schnablelab import __version__ as version
from schnablelab.apps.Tools import print_json
from schnablelab.apps.base import OptionParser, OptionGroup, ActionDispatcher, SUPPRESS_HELP
from multiprocessing.dummy import Pool as ThreadPool
from retrying import retry

# API key and URLs
PLANET_API_KEY = os.getenv('PLANET_API_KEY')
URL = "https://api.planet.com/data/v1"

default_geojson = op.abspath(op.dirname(__file__)) + '/lindsay_james.geojson'

session = requests.Session()
session.auth = (PLANET_API_KEY, '')

def mapdict2df(dict_list):
    '''
    loop dicts in a list and return a dataframe with each row representing a dict
    '''
    tmp_df = []
    for i in dict_list:
        tmp_df.append(pd.DataFrame(i))
    df = pd.concat(tmp_df)
    return df

class Client(object):
    """
    wrapper of planet API
    """
    def __init__(self, url=URL, ses=session):
        self.url = url
        self.ses = session
        self.url_stat = "{}/stats".format(self.url)
        self.url_search = "{}/searches".format(self.url)
        self.url_quick_search = "{}/quick-search".format(self.url)
        
    def get_all_items(self):
        res = self.ses.get(self.url)
        if res.status_code == 200:
            url_item_types = res.json()['_links']['item-types']
            item_types_list = self.ses.get(url_item_types).json()['item_types']
            item_types_df = mapdict2df(item_types_list)
            return item_types_df

    def get_all_assets(self):
        res = self.ses.get(self.url)
        if res.status_code == 200:
            url_asset_types = res.json()['_links']['asset-types']
            item_asset_list = self.ses.get(url_asset_types).json()['asset_types']
            item_asset_df = mapdict2df(item_asset_list)
            return item_asset_df

def parse_date(str_y_m_d):
    y, m, d = str_y_m_d.split('-')
    return "%s-%s-%sT00:00:00.000Z"%(y, m, d)

class Filter(object):
    '''
    customize filter object
    '''
    def __init__(self):
        pass

    def get_filter_date(self, field_name='acquired', st='', ed=''):
        date_dict, config_dict = dict(), dict()
        if st:
            config_dict['gte'] = parse_date(st)
        if ed:
            config_dict['lte'] = parse_date(ed)
        date_dict['type'] = 'DateRangeFilter'
        date_dict['field_name'] = field_name
        date_dict['config'] = config_dict
        return date_dict
    
    def get_filter_cloud(self, st=0, ed=0.05):
        return {
            'type': 'RangeFilter',
            'field_name': 'cloud_cover',
            'config': {'gte': st, 'lte': ed}
        }
    
    def get_filter_coverage(self, st=0, ed=1):
        return {
            'type': 'RangeFilter',
            'field_name': 'usable_data',
            'config': {'gte': st, 'lte': ed}
        }

    def get_filter_instrument(self, instrument=['PS2']):
        return {
            "type": "StringInFilter",
            "field_name": "instrument",
            "config": instrument
        }
    
    def get_filter_permission(self):
        return {
            "type": "PermissionFilter",
            "config": ["assets.analytic:download"]
        }

    def get_filter_geometry(self, config):
        return {
            'type': 'GeometryFilter',
            'field_name': 'geometry',
            'config': config
        }

class Request(object):
    '''
    customzie request object
    stat endpoint requires the interval
    search endpoint doesn't require the interval
    '''
    def __init__(self, item_types = ['PSScene4Band', 'PSOrthoTile'], interval=''):
        self.item_types = item_types
        self.interval = interval if interval else None
    
    def get_request(self, filter):
        request_dict = dict()
        request_dict['item_types'] = self.item_types
        if self.interval:
            request_dict['interval'] = self.interval
        request_dict['filter'] = filter
        return request_dict
        
def asset_types(args):
    '''
    %prog asset_types 

    print all available asset types
    '''
    p = OptionParser(asset_types.__doc__)
    p.add_option('-o', '--output', default="asset_types.csv",
                help='specify output file')
    p.add_option('--n', default="all",
                help='how many rows you wanna see')
    opts, args = p.parse_args(args)
    if len(args) != 0:
        sys.exit(not p.print_help())
    
    client = Client()
    df_assets = client.get_all_assets()[['id', 'display_name', 'display_description']]
    if opts.n != 'all':
        try:
            rows = int(opts.n)
            df_assets.head(rows).to_csv(opts.output, index=False, sep='\t')
        except ValueError:
            sys.exit("--n must be a number")
    else:
        df_assets.to_csv(opts.output, index=False, sep='\t')
    print('check %s!'%(opts.output))


def item_types(args):
    '''
    %prog item_types
    print all available item types
    '''
    p = OptionParser(item_types.__doc__)
    p.add_option('-o', '--output', default="item_types.csv",
                help='specify output file')
    p.add_option('--n', default="all",
                help='how many rows you wanna see')
    opts, args = p.parse_args(args)
    if len(args) != 0:
        sys.exit(not p.print_help())
    
    client = Client()
    df_items = client.get_all_items()[['id', 'display_name', 'display_description']]
    if opts.n != 'all':
        try:
            rows = int(opts.n)
            df_items.head(rows).to_csv(opts.output, index=False, sep='\t')
        except ValueError:
            sys.exit("n must be a number")
    else:
        df_items.to_csv(opts.output, index=False, sep='\t')
    print('check %s!'%(opts.output))

def stat(args):
    '''
    %prog stat
    
    Check availabe images on Planet
    '''
    p = OptionParser(stat.__doc__)
    p.add_option('-o', '--output', default="stats.csv",
                help='specify output file')
    p.add_option('--geom', default="lindsay_james.geojson",
                help='speficy the geojson file containing the geometry info')
    p.add_option('--cloud', default=True, action='store_false',
                help = 'disable cloud filter if add --cloud option')
    p.add_option('--coverage', default=False, action='store_true',
                help = 'add area coverage filter if add --coverage option')
    p.add_option('--instrument', default=False, action='store_true',
                help = 'add instrument filter if add --instrument option')
    p.add_option('--date_range', default=True, action='store_false',
                help = 'disable date filter if add --date_range option')
    
    q = OptionGroup(p, "options for date filter")
    p.add_option_group(q)
    q.add_option('--start', default="2018-01-01",
                help='the start date. use yyyy-mm-dd format.')
    q.add_option('--end',
                help='the end date. use yyyy-mm-dd format')
    
    r = OptionGroup(p, "options for request")
    p.add_option_group(r)
    r.add_option('--interval', default="year", choices=('year', 'month', 'week', 'day', 'hour'),
                help='specify the interval in the request')
    r.add_option('--item_types', default='PSScene4Band,PSOrthoTile', 
                help='specify the item types. use comma separated if more than one item')

    opts, args = p.parse_args(args)
    if len(args) != 0:
        sys.exit(not p.print_help())

    # all filters
    Filters, Filters_names = [], []
    geojson_fn = default_geojson if opts.geom == 'lindsay_james.geojson' else opts.geom
    with open(geojson_fn) as f:
        data = json.load(f)
    geometry = data['features'][0]['geometry']
    filter_geom = Filter().get_filter_geometry(geometry)
    Filters.append(filter_geom)
    Filters_names.append(filter_geom['field_name'])
    if opts.cloud:
        filter_cloud = Filter().get_filter_cloud()
        Filters.append(filter_cloud)
        Filters_names.append(filter_cloud['field_name'])
    if opts.coverage:
        filter_coverage = Filter().get_filter_coverage()
        Filters.append(filter_coverage)
        Filters_names.append(filter_coverage['field_name'])
    if opts.instrument:
        filter_instrument = Filter().get_filter_instrument()
        Filters.append(filter_instrument)
        Filters_names.append(filter_instrument['field_name'])
    if opts.date_range:
        st, ed = opts.start, opts.end
        filter_date = Filter().get_filter_date(st=st, ed=ed)
        Filters.append(filter_date)
        Filters_names.append(filter_date['field_name'])
    Final_Filters = {
        'type': 'AndFilter',
        'config': Filters
    }
    print('Applied Filters: %s'%(', '.join(Filters_names)))
    
    # request
    Requests = Request(item_types = opts.item_types.split(','), interval=opts.interval)
    Final_Requests = Requests.get_request(Final_Filters)    

    # post
    client = Client()
    res = client.ses.post(client.url_stat, json=Final_Requests)
    df_stat = pd.DataFrame.from_dict(res.json()['buckets'])[['start_time', 'count']]
    print(df_stat)
    df_stat.to_csv(opts.output, index=False, sep='\t')
    print('also saved to %s!'%opts.output)

def quick_search(args):
    '''
    %prog quick_search
    
    Perform quick serach and get ids for items from search resutls
    '''
    p = OptionParser(quick_search.__doc__)
    p.add_option('-o', '--output', default="searches.csv",
                help='specify output file')
    p.add_option('--geom', default="lindsay_james.geojson",
                help='speficy the geojson file containing the geometry info')
    p.add_option('--cloud', default=True, action='store_false',
                help = 'disable cloud filter if add --cloud option')
    p.add_option('--coverage', default=False, action='store_true',
                help = 'add area coverage filter if add --coverage option')
    p.add_option('--instrument', default=False, action='store_true',
                help = 'add instrument filter if add --instrument option')
    p.add_option('--date_range', default=True, action='store_false',
                help = 'disable date filter if add --date_range option')
    p.add_option('--map_footprint', default=False, action='store_true',
                help = 'add mapping footprints if add --map_footprint option')
    
    q = OptionGroup(p, "options for date filter")
    p.add_option_group(q)
    q.add_option('--start', default="2018-01-01",
                help='the start date. use yyyy-mm-dd format.')
    q.add_option('--end',
                help='the end date. use yyyy-mm-dd format')
    
    r = OptionGroup(p, "options for request")
    p.add_option_group(r)
    r.add_option('--item_types', default='PSScene4Band', 
                help='specify the item types. use comma separated if more than one item')

    opts, args = p.parse_args(args)
    if len(args) != 0:
        sys.exit(not p.print_help())

    # all filters
    Filters, Filters_names = [], []
    geojson_fn = default_geojson if opts.geom == 'lindsay_james.geojson' else opts.geom
    with open(geojson_fn) as f:
        data = json.load(f)
    geometry = data['features'][0]['geometry']
    filter_geom = Filter().get_filter_geometry(geometry)
    Filters.append(filter_geom)
    Filters_names.append(filter_geom['field_name'])
    if opts.cloud:
        filter_cloud = Filter().get_filter_cloud()
        Filters.append(filter_cloud)
        Filters_names.append(filter_cloud['field_name'])
    if opts.coverage:
        filter_coverage = Filter().get_filter_coverage()
        Filters.append(filter_coverage)
        Filters_names.append(filter_coverage['field_name'])
    if opts.instrument:
        filter_instrument = Filter().get_filter_instrument()
        Filters.append(filter_instrument)
        Filters_names.append(filter_instrument['field_name'])
    if opts.date_range:
        st, ed = opts.start, opts.end
        filter_date = Filter().get_filter_date(st=st, ed=ed)
        Filters.append(filter_date)
        Filters_names.append('date_range')
    Final_Filters = {
        'type': 'AndFilter',
        'config': Filters
    }
    print('Applied Filters: %s'%(', '.join(Filters_names)))
    
    # request
    items = opts.item_types.split(',')
    print('Item Types: ', items)
    Requests = Request(item_types = items)
    Final_Requests = Requests.get_request(Final_Filters)    
    #rint(Final_Requests)

    # post
    client = Client()
    res = client.ses.post(client.url_quick_search, json=Final_Requests, params={"_page_size" : 10})
    geojson = res.json() # important keys: _links(current link and link of the next page), features(id)
    link_first_page = geojson['_links']['_first']

    ids, cloud_cover, item_type, assets_url = [], [], [], []
    def parse_page(session, search_url, map_footprint):
        '''
        loop pages and extract ids from each page
        '''
        res = session.get(search_url)
        if map_footprint:
            url = geojsonio.display(res.text)
        page = res.json()
        for feature in page['features']:
            id = feature['id']
            cc = feature['properties']['cloud_cover']
            it = feature['properties']['item_type']
            au = feature['_links']['assets']
            ids.append(id)
            cloud_cover.append(cc)
            item_type.append(it)
            assets_url.append(au)
        next_url = page["_links"].get("_next")
        if next_url:
            parse_page(session, next_url, map_footprint)
    parse_page(client.ses, link_first_page, opts.map_footprint)

    df = pd.DataFrame(dict(zip(['id', 'cloud_cover', 'item_type', 'assets_url'], [ids, cloud_cover, item_type, assets_url])))
    df.to_csv(opts.output, index=False, sep='\t')
    print('%s items found, please check search resutls in %s!'%(df.shape[0], opts.output))

def activate(args):
    '''
    %prog activate searches.csv

    each id represt an item, each item has many different assets listed in the assets endpoint, 
    such as analytic and analytic_xml...
    This function will perform activation for the selecetd asset so they can be downloaded later 
    '''
    p = OptionParser(activate.__doc__)
    p.add_option('--output', default="download_links_AssetKey.csv",
                help='specify output file')
    p.add_option('--asset_key', default="analytic_xml", choices=('analytic', 'analytic_xml', 'visual'),
                help='speficy the asset you are goting to activate and download')
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())
    searches, = args

    client = Client()
    asset_urls_df = pd.read_csv(searches, delim_whitespace=True)
    asset_urls_df['asset_type'] = opts.asset_key

    
    # perform activation (control request rate here)
    codes = []
    for idx, row in asset_urls_df.iterrows():
        print(idx, row['id'], row['item_type'])
        res = client.ses.get(row['assets_url']) # normal request
        activate_url = res.json()[opts.asset_key]['_links']['activate']
        res = client.ses.get(activate_url) # activate request
        code = res.status_code
        codes.append(code)
        if code == 202:
            print('activation started')
        elif code == 204:
            print('activated')
        elif code == 401:
            print('no permission')
        else:
            print(code)
        time.sleep(.500)
    asset_urls_df['activation_code'] = np.array(codes)
    
    # check activation status (control request rate here)
    asset_activated = False
    while asset_activated == False:
        if 202 not in asset_urls_df['activation_code'].values:
            asset_activated = True
            print('all assets have been activated!!!')
        else:
            asset_urls_df_202 = asset_urls_df.copy()
            for idx, row in asset_urls_df_202.iterrows():
                if row['activation_code'] == 202:
                    status = client.ses.get(row['assets_url']).json()[opts.asset_key]['status'] # normal request
                    time.sleep(.200)
                    if status == 'active':
                        print(row['id'], 'activated!')
                        asset_urls_df.loc[idx, 'activation_code'] = 204
                    else: print(row['id'], 'activating...')
        time.sleep(5) # sleep 5s to start another round check

    # update link for downloading
    def get_download_link(asset_link):
        res = client.ses.get(asset_link).json()[opts.asset_key] # normal request
        time.sleep(.200)
        return res['location'] if res['status'] == 'active' else np.nan
    asset_urls_df['download_link'] = asset_urls_df['assets_url'].apply(get_download_link)

    output = opts.output.replace('AssetKey', opts.asset_key) \
        if opts.output == "download_links_AssetKey.csv" \
        else opts.output
    asset_urls_df.to_csv(output, index=False, sep='\t', na_rep='NaN')
    print('check the download links in %s!'%output)

def download(args):
    '''
    %prog activate download_links.csv

    download activated asset links 
    '''
    p = OptionParser(download.__doc__)
    p.add_option('--output', default="'infer'",
                help='default to construct the output file name from the API response')
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(not p.print_help())
    links_csv, = args
    links_df = pd.read_csv(links_csv, delim_whitespace=True)
    client = Client()

    # Send a GET request to the provided location url,
    for __, row in links_df.iterrows():
        res = client.ses.get(row['download_link'], stream=True)
        suffix = 'tif' if row['asset_type']=='visual' else row['asset_type']
        output = '%s_%s.%s'%(row['id'], row['item_type'], suffix) \
            if opts.output == "'infer'" else opts.output
        # Save the file
        with open(output, "wb") as f:
            print('download %s...'%output)
            for chunk in res.iter_content(chunk_size=1024):
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    f.flush()


def main():
    actions = (
        ('asset_types', 'print all available asset types'),
        ('item_types', 'print all available item types'),
        ('stat', 'check available imagery counts'),
        ('quick_search', 'perform quick search to get all the target ids'),
        ('activate', 'activate assets for downloading'),
        ('download', 'download activated links')
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

if __name__ == '__main__':
    main()
