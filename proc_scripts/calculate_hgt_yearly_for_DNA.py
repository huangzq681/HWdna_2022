'''
This script is for calculating JJA geopotential height anomalies. GPH anomalies for reanalyses and CMIP6 forcing are computed by 
removing the multi-year average from the yearly GPH series.
'''
import xarray as xr
import pandas as pd
import numpy as np
import os
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

dataset_src_run = {
    'CanESM5':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'piControl':['r1i1p1f1']
    },
    'HadGEM3-GC31-LL':{
        'historical':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-GHG':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-nat':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-aer':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'ssp585':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'piControl':['r1i1p1f1']
    },
    'MIROC6':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'piControl':['r1i1p1f1']
    },
    'IPSL-CM6A-LR':{
        'historical':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r6i1p1f1'], #'r3i1p1f1','r14i1p1f1' are excluded, as there are not corresponding runs in historical forcing
        'piControl':['r1i2p1f1']
    },
    'MRI-ESM2-0':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'piControl':['r1i1p1f1']
    },
}


time_range_hgt = {
    'era5':'1950-2020',
    'jra55':'1958-2014',
    'ncep2':'1979-2020',
    'CanESM5':{
        'historical':'19510101-20141231',
        'hist-GHG':'19510101-20201231',
        'hist-nat':'19510101-20201231',
        'hist-aer':'19510101-20201231',
        'ssp585':'20510101-21001231',
        'piControl':'52010101-62001231'
    },
    'HadGEM3-GC31-LL':{
        'historical':'19500101-20141230',
        'hist-GHG':'19500101-20201230',
        'hist-nat':'19500101-20201230',
        'hist-aer':'19500101-20201230',
        'ssp585':'20500101-20991230',
        'piControl':'18500101-23491230'
    },
    'MIROC6':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20500101-21001231',
        'piControl':'32000101-36991231'
    },
    'IPSL-CM6A-LR':{
        'historical':'19500101-20141231',
        'hist-GHG':'19600101-20201231',
        'hist-nat':'19600101-20201231',
        'hist-aer':'19600101-20201231',
        'ssp585':'20150101-21001231',
        'piControl':'18700101-20991231'
    },
    'MRI-ESM2-0':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20550101-21001231',
        'piControl':'18500101-20491231'
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
    var_name = 'hgt'
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    var_name_tran = {'hgt':'tasmax','hgt':'zg'}
    if forcing in ['hist-GHG','hist-aer','hist-nat','ssp585','piControl']:
        if dataset_name != 'IPSL-CM6A-LR':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
        else:
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
    else:
        if dataset_name in ['era5','jra55','ncep2']:
            filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_hgt[dataset_name] + '.nc'
        else:
            if dataset_name != 'IPSL-CM6A-LR':
                filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
            else:
                filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
    if dataset_name in ['era5','jra55','ncep2']:
        data = xr.open_dataarray(filepath)
        if 'height' in data.dims:
            data = data.squeeze(drop=True)
    else:
        data = xr.open_dataset(filepath)
        data = data[var_name_tran[var_name]]
    if dataset_name == 'era5':
        data = data.rename({'longitude':'lon','latitude':'lat','time':'time'})
    if forcing != 'ssp585' and forcing != 'piControl': 
        if dataset_name == 'HadGEM3-GC31-LL':
            time_start = '1979-01-01'
            time_end = '2014-12-30'
        else:
            time_start = '1979-01-01'
            time_end = '2014-12-31'
    elif forcing == 'ssp585':
        if dataset_name == 'HadGEM3-GC31-LL':
            time_start = '2064-01-01'
            time_end = '2099-12-30'
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

## Calculate GPH anomalies
def cal_gph_anomalies(dataset_name,domain,forcing=None,run=None):
    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    if dataset_name in ['era5','jra55','ncep2']:
        rd_dat = sel_domain(dataset_name=dataset_name,domain=domain)
    else:
        rd_dat = sel_domain(dataset_name=dataset_name,domain=domain,forcing=forcing,run=run)
    
    weights = np.cos(np.deg2rad(rd_dat.lat))
    weights.name = "weights"

    hgt_occur_weight = rd_dat.weighted(weights)
    hgt_occur_weight_mean = hgt_occur_weight.mean(("lon", "lat"))
    hgt_occur_weight_mean_yearly = hgt_occur_weight_mean.groupby('time.year').mean()
    hgt_occur_weight_mean_yearly = hgt_occur_weight_mean_yearly.squeeze(drop=True)

    # hgt_summer_avg = hgt_occur_weight_mean_yearly.mean(dim='time')
    hgt_ds_summer = hgt_occur_weight_mean_yearly #- hgt_summer_avg
    
    return hgt_ds_summer


def save_result(domain):
    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    hgt_occur_reanalyses = pd.DataFrame(columns=['era5','jra55','ncep2'])

    path = to_file_dir + 'GPH_reanlyses_yearly_all_patterns_' + domain + '.csv'
    if os.path.exists(path) == True:
        pass
    else:
        for d in ['era5','jra55','ncep2']:
            hgt_occur = cal_gph_anomalies(dataset_name=d,domain=domain)
            hgt_occur = hgt_occur.to_series()
            print(hgt_occur)
            # hgt_occur.index = range(1979,2015)
            hgt_occur_reanalyses[d] = hgt_occur
            print(d)
        hgt_occur_reanalyses.to_csv(to_file_dir + 'GPH_reanlyses_yearly_all_patterns_' + domain + '.csv')


    for f in ['historical','hist-GHG','hist-aer','hist-nat']:
        path = to_file_dir + 'GPH_' + f + '_yearly_all_patterns_' + domain + '.csv'
        if os.path.exists(path) == True:
            pass 
        else:
            data_run_name = []
            for d in dataset_src_run.keys():
                for k in dataset_src_run[d][f]:
                    data_run_name.append(d + '_' + k)
            hgt_occur_cmip = pd.DataFrame(columns=data_run_name)
            for d in dataset_src_run.keys():
                for r in dataset_src_run[d][f]:
                    hgt_occur = cal_gph_anomalies(dataset_name=d,domain=domain,forcing=f,run=r)
                    hgt_occur = hgt_occur.to_series()
                    print(hgt_occur)
                    # hgt_occur.index = range(1979,2015)
                    hgt_occur_cmip[d + '_' + r] = hgt_occur
                    print(f + '_' + d + '_' + r)
            hgt_occur_cmip.to_csv(path)
    
save_result('WNA')
save_result('EU')
save_result('EAS')

