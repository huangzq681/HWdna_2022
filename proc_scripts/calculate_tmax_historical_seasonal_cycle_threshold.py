'''
This script is for calculating JJA hot extreme occurrence. Temperature anomalies for reanalyses and CMIP6 historical forcing are computed by 
removing the seasonal cycle from daily reanalysis 2-m maximum temperature. Temperature anomalies for CMIP6 external forcings are computed by
removing the historical seasonal cycle from daily reanalysis 2-m maximum temperature of the same model run.
Hot extreme thresholds are defined as the 95th percentile value of the 1979-2014 daily 2-m maximum temperature anomaly distribution.
Hot/cold extreme occurrences are defined as days on which the daily temperature anomalies are greater (or equal to) the hot extreme thresholds.
'''
import xarray as xr
import pandas as pd
import numpy as np

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
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r6i1p1f1','r14i1p1f1'],
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
    'era5':'1950-2020',
    'jra55':'1958-2014',
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
    # 'WNA':{'lon_min':230,'lon_max':262,'lat_min':30,'lat_max':54},
    'WNA':{'lon_min':220,'lon_max':260,'lat_min':25,'lat_max':55},
}

def sel_lonlat_range(dataarray,lon_min,lon_max,lat_min,lat_max):
    mask_lon = (dataarray.lon >= lon_min) & (dataarray.lon <= lon_max)
    mask_lat = (dataarray.lat >= lat_min) & (dataarray.lat <= lat_max)
    dataarray = dataarray.where(mask_lon & mask_lat, drop=True)
    return dataarray

def sel_domain(dataarray,lon_min,lon_max,lat_min,lat_max):
    mask_lon = (dataarray.lon >= lon_min) & (dataarray.lon <= lon_max)
    mask_lat = (dataarray.lat >= lat_min) & (dataarray.lat <= lat_max)
    dataarray = dataarray.where(mask_lon & mask_lat, drop=True)
    return dataarray

def season_cycle_tmax_threshold(dataset_name,domain,forcing=None,run=None):
    var_name = 'tmax'
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    var_name_tran = {'hgt':'zg','tmax':'tasmax'}
    if dataset_name in ['era5','jra55','ncep2']:
        filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_tmax[dataset_name] + '.nc'
    else:
        if dataset_name != 'IPSL-CM6A-LR':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_tmax[dataset_name][forcing] + '.nc'
        else:
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_tmax[dataset_name][forcing] + '.nc'
    if dataset_name in ['era5','jra55','ncep2']:
        data = xr.open_dataarray(filepath)
    else:
        data = xr.open_dataset(filepath)
        data = data[var_name_tran[var_name]]
    if dataset_name == 'era5':
        data = data.rename({'longitude':'lon','latitude':'lat','time':'time'})
        
    if dataset_name == 'HadGEM3-GC31-LL':
        time_start = '1979-01-01'
        time_end = '2014-12-30'
    else:
        time_start = '1979-01-01'
        time_end = '2014-12-31'
        
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

    tmax_season_cycle_quantile = data_summer.groupby('time.dayofyear').quantile(q=0.95)

    if dataset_name in ['era5','jra55','ncep2']:
        to_file_name = '/home/xtan/scratch/hzq/HWdna/procData/'+ dataset_name + '_' + var_name + '_season_cycle_quantile_1979-2014' + domain + '.nc'
    else:
        to_file_name = '/home/xtan/scratch/hzq/HWdna/procData/'+ dataset_name + '_' + forcing + '_' + run + '_' + var_name + '_season_cycle_quantile_1979-2014' + domain + '.nc'
    tmax_season_cycle_quantile.to_netcdf(to_file_name)

def save_result(domain):
    for d in time_range_tmax.keys():
        if d in ['era5','jra55','ncep2']:
            season_cycle_tmax_threshold(d,domain=domain)
            print(d)
        else:
            for f in ['historical']:
                for r in dataset_src_run[d][f]:
                    season_cycle_tmax_threshold(dataset_name=d,domain=domain,forcing=f,run=r)
                    print(d+'-'+f+'-'+r)
save_result('WNA')
save_result('EU')
save_result('EAS')