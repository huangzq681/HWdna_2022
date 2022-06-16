'''
This script is for selecting 500hpa geopotential height corresponding to each patterns and calculating their averages over the research periods
from CMIIP6 models under different external forcings
'''

import xarray as xr
import pandas as pd
import numpy as np
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

dataset_src_run = {
    'CanESM5':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
    },
    'HadGEM3-GC31-LL':{
        'historical':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-GHG':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-nat':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-aer':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'ssp585':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
    },
    'MIROC6':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
    },
    'IPSL-CM6A-LR':{
        'historical':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r6i1p1f1'], #'r3i1p1f1','r14i1p1f1' are excluded, as there are not corresponding runs in historical forcing
    },
    'MRI-ESM2-0':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
    },
}

time_range_tmax = {
    'jra55':'1958-2014',
    'era5':'1950-2020',
    'ncep2':'1979-2020',
    'CanESM5':{
        'historical':'18500101-20141231',
        'hist-GHG':'18500101-20201231',
        'hist-nat':'18500101-20201231',
        'hist-aer':'18500101-20201231',
        'ssp585':'20150101-21001231'
    },
    'HadGEM3-GC31-LL':{
        'historical':'19500101-20141230',
        'hist-GHG':'19500101-20201230',
        'hist-nat':'19500101-20201230',
        'hist-aer':'19500101-20201230',
        'ssp585':'20500101-21001230'
    },
    'MIROC6':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20191231',
        'hist-nat':'19500101-20191231',
        'hist-aer':'19500101-20191231',
        'ssp585':'20550101-21001231'
    },
    'IPSL-CM6A-LR':{
        'historical':'18500101-20141231',
        'hist-GHG':'18500101-20201231',
        'hist-nat':'18500101-20201231',
        'hist-aer':'18500101-20201231',
        'ssp585':'20150101-21001231'
    },
    'MRI-ESM2-0':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20650101-21001231'
    }
}

domain_lonlat = {
    'EAS':{'lon_min':90,'lon_max':130,'lat_min':30,'lat_max':60},
    'EU':{'lon_min':10,'lon_max':50,'lat_min':35,'lat_max':65},
    'WNA':{'lon_min':220,'lon_max':260,'lat_min':25,'lat_max':55},
}

target_griddes_sub = {
    'EAS':{
        'lat': np.arange(domain_lonlat['EAS']['lat_min'], domain_lonlat['EAS']['lat_max']+1, 1),
        'lon':np.arange(domain_lonlat['EAS']['lon_min'], domain_lonlat['EAS']['lon_max']+1, 1)},
    'EU':{
        'lat': np.arange(domain_lonlat['EU']['lat_min'], domain_lonlat['EU']['lat_max']+1, 1),
        'lon':np.arange(domain_lonlat['EU']['lon_min'], domain_lonlat['EU']['lon_max']+1, 1)},
    'WNA':{
        'lat': np.arange(domain_lonlat['WNA']['lat_min'], domain_lonlat['WNA']['lat_max']+1, 1),
        'lon':np.arange(domain_lonlat['WNA']['lon_min'], domain_lonlat['WNA']['lon_max']+1, 1)}
}

target_griddes = {'lat': np.arange(0, 90, 1.5),'lon':np.arange(0, 360, 1.5)}

def sel_lonlat_range(dataarray,lon_min,lon_max,lat_min,lat_max):
    mask_lon = (dataarray.lon >= lon_min) & (dataarray.lon <= lon_max)
    mask_lat = (dataarray.lat >= lat_min) & (dataarray.lat <= lat_max)
    dataarray = dataarray.where(mask_lon & mask_lat, drop=True)
    return dataarray

def sel_domain(dataset_name,domain,forcing=None,run=None):
    var_name = 'tmax'
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    var_name_tran = {'tmax':'tasmax','hgt':'zg'}
    if forcing in ['hist-GHG','hist-aer','hist-nat','ssp585']:
        if dataset_name != 'IPSL-CM6A-LR':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_tmax[dataset_name][forcing] + '.nc'
        else:
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_tmax[dataset_name][forcing] + '.nc'
    else:
        if dataset_name in ['era5','jra55','ncep2']:
            filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_tmax[dataset_name] + '.nc'
        else:
            if dataset_name != 'IPSL-CM6A-LR':
                filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_tmax[dataset_name][forcing] + '.nc'
            else:
                filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_tmax[dataset_name][forcing] + '.nc'
    if dataset_name in ['era5','jra55','ncep2']:
        data = xr.open_dataarray(filepath)
        if 'height' in data.dims:
            data = data.squeeze(drop=True)
    else:
        data = xr.open_dataset(filepath)
        data = data[var_name_tran[var_name]]
    if dataset_name == 'era5':
        data = data.rename({'longitude':'lon','latitude':'lat','time':'time'})
    if forcing != 'ssp585': 
        if dataset_name == 'HadGEM3-GC31-LL':
            time_start = '1979-01-01'
            time_end = '2014-12-30'
        else:
            time_start = '1979-01-01'
            time_end = '2014-12-31'
    else:
        if dataset_name == 'HadGEM3-GC31-LL':
            time_start = '2064-01-01'
            time_end = '2099-12-30'
        elif dataset_name == 'MRI-ESM2-0':
            time_start = '2065-01-01'
            time_end = '2099-12-31'
        else:
            time_start = '2064-01-01'
            time_end = '2099-12-31'
    data = data.sel(time=slice(time_start,time_end))
    data_summer = data.sel(time=data['time.season']=='JJA')
    if domain == 'EAS':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['EAS']['lon_min'],lon_max=domain_lonlat['EAS']['lon_max'],lat_min=domain_lonlat['EAS']['lat_min'],lat_max=domain_lonlat['EAS']['lat_max'])
    elif domain == 'EU':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['EU']['lon_min'],lon_max=domain_lonlat['EU']['lon_max'],lat_min=domain_lonlat['EU']['lat_min'],lat_max=domain_lonlat['EU']['lat_max'])
    elif domain == 'WNA':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['WNA']['lon_min'],lon_max=domain_lonlat['WNA']['lon_max'],lat_min=domain_lonlat['WNA']['lat_min'],lat_max=domain_lonlat['WNA']['lat_max'])
    else:
        raise Exception('ERROR, domain must be a string object belongs to [\'NH\',\'EAS\',\'EU\',\'WNA\']')
    return data_summer

def to_ncfile(domain,forcing):
    to_dir = '/lustre07/scratch/xtan/hzq/HWdna/procData/' + 'tmax'+ '_average_' + forcing + '_' + domain + '.nc'
    return to_dir

def calculate_tmax_avg(forcing,domain):

    dataset_name = 'CanESM5'
    run = 'r1i1p1f1'
    rd_dat = sel_domain(dataset_name=dataset_name,domain=domain,forcing=forcing,run=run)
    rd_dat = rd_dat.mean('time')

    rd_dat = rd_dat.interp(lat=target_griddes_sub[domain]['lat'],lon=target_griddes_sub[domain]['lon'])
    rd_dat = rd_dat.interpolate_na(dim='lon',method='linear',fill_value='extrapolate')
    rd_dat = rd_dat.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
    rd_dat_all = rd_dat.expand_dims({'Dataset_run':1})
    rd_dat_all = rd_dat_all.assign_coords(Dataset_run = ['CanESM5_r1i1p1f1'])
    
    for dataset_name in dataset_src_run.keys():
        for run in dataset_src_run[dataset_name][forcing]:
            if dataset_name == 'CanESM5' and run == 'r1i1p1f1':
                pass
            else:
                rd_dat = sel_domain(dataset_name=dataset_name,domain=domain,forcing=forcing,run=run)
                rd_dat = rd_dat.mean('time')
                rd_dat = rd_dat.interp(lat=target_griddes_sub[domain]['lat'],lon=target_griddes_sub[domain]['lon'])
                rd_dat = rd_dat.interpolate_na(dim='lon',method='linear',fill_value='extrapolate')
                rd_dat = rd_dat.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
                rd_dat = rd_dat.expand_dims({'Dataset_run':1})
                rd_dat = rd_dat.assign_coords(Dataset_run = [dataset_name + '_' + run])
                rd_dat_all = xr.concat([rd_dat_all,rd_dat],dim='Dataset_run')

    rd_dat_all.to_netcdf(to_ncfile(domain=domain,forcing=forcing))

for f in ['historical']:
    for d in ['EU','EAS','WNA']:
        calculate_tmax_avg(forcing=f,domain=d)
        print(f + '_' + d)

        