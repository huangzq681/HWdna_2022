'''
This script is for calculating JJA hot extreme occurrence for piControl. Temperature anomalies for reanalyses and CMIP6 historical forcing are computed by 
removing the seasonal cycle from daily reanalysis 2-m maximum temperature. Temperature anomalies for CMIP6 external forcings are computed by
removing the historical seasonal cycle from daily reanalysis 2-m maximum temperature of the same model run.
Hot extreme thresholds are defined as the 95th percentile value of the 1979-2014 daily 2-m maximum temperature anomaly distribution.
Hot/cold extreme occurrences are defined as days on which the daily temperature anomalies are greater (or equal to) the hot extreme thresholds.
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
        'piControl':['r1i1p1f3']
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
        'piControl':['r1i1p1f1']
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

MIROC6_ts = [str(3200 + i * 10) + '0101-' + str(3200 + i * 10 + 9 ) + '1231' for i in range(50)]

time_range_tmax = {
    'jra55':'1958-2014',
    'era5':'1950-2020',
    'ncep2':'1979-2020',
    'CanESM5':{
        'historical':'18500101-20141231',
        'hist-GHG':'18500101-20201231',
        'hist-nat':'18500101-20201231',
        'hist-aer':'18500101-20201231',
        'ssp585':'20150101-21001231',
        # 'piControl':'52000101-62001231'
        'piControl':['52010101-54001231','54010101-56001231','56010101-58001231','58010101-60001231','60010101-62001231']
    },
    'HadGEM3-GC31-LL':{
        'historical':'19500101-20141230',
        'hist-GHG':'19500101-20201230',
        'hist-nat':'19500101-20201230',
        'hist-aer':'19500101-20201230',
        'ssp585':'20500101-21001230',
        # 'piControl':'18500101-20491230'
        'piControl':['18500101-19491230','19500101-20491230','20500101-21491230','21500101-22491230','22500101-23491230']
    },
    'MIROC6':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20191231',
        'hist-nat':'19500101-20191231',
        'hist-aer':'19500101-20191231',
        'ssp585':'20550101-21001231',
        # 'piControl':'32000101-36991231'
        'piControl': MIROC6_ts
    },
    'IPSL-CM6A-LR':{
        'historical':'18500101-20141231',
        'hist-GHG':'18500101-20201231',
        'hist-nat':'18500101-20201231',
        'hist-aer':'18500101-20201231',
        'ssp585':'20150101-21001231',
        # 'piControl':'18500101-38491231',
        'piControl':['18500101-23491231','23500101-28491231','28500101-30491231','30500101-35491231','35500101-38491231'],
    },
    'MRI-ESM2-0':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20650101-21001231',
        # 'piControl':'18500101-20491231',
        'piControl':['18500101-18991231','19000101-19491231','19500101-19991231','20000101-20491231']
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

def sel_domain(dataset_name,domain,forcing=None,run=None,time=None):
    var_name = 'tmax'
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    var_name_tran = {'tmax':'tasmax'}
    if dataset_name != 'IPSL-CM6A-LR':
        filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time + '.nc'
    else:
        filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time + '.nc'
    data = xr.open_dataset(filepath)
    data = data[var_name_tran[var_name]]
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

## Calculate Hot extreme occurrences as days on which the daily temperature anomalies are greater (or equal to) the hot extreme thresholds
## season_cycle_data is the seasonal cycle from daily reanalysis 2-m maximum temperature
def cal_hot_extreme_anomalies(regrid_domain_dat,season_cycle_data,domain):
    tmax_ds_summer = regrid_domain_dat
    tmax_ds_summer = tmax_ds_summer.groupby('time.dayofyear') - season_cycle_data
    tmax_ds_summer_ishot = (tmax_ds_summer > 0) * 1
    tmax_ds_summer = tmax_ds_summer * tmax_ds_summer_ishot
    # regrid
    tmax_ds_summer = tmax_ds_summer.interp(lat=target_griddes_sub[domain]['lat'],lon=target_griddes_sub[domain]['lon'])
    tmax_ds_summer = tmax_ds_summer.interpolate_na(dim='lon',method='linear',fill_value='extrapolate')
    tmax_ds_summer = tmax_ds_summer.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
    return tmax_ds_summer


def cal_hot_extreme_occur(dataset_name,domain,forcing=None,run=None,time=None):
    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    rd_dat = sel_domain(dataset_name=dataset_name,domain=domain,forcing=forcing,run=run,time=time)
    season_cycle_quantile_path = to_file_dir + dataset_name + '_historical_' + run + '_tmax_season_cycle_quantile_1979-2014' + domain + '.nc'
    season_cycle_quantile = xr.open_dataarray(season_cycle_quantile_path)
    tmax_ano = rd_dat.groupby('time.dayofyear') - season_cycle_quantile
    tmax_occur = tmax_ano > 0 * 1
    weights = np.cos(np.deg2rad(tmax_ano.lat))
    weights.name = "weights"
    tmax_occur_weight = tmax_occur.weighted(weights)
    tmax_occur_weight_mean = tmax_occur_weight.mean(("lon", "lat"))
    tmax_occur_weight_mean_yearly = tmax_occur_weight_mean.groupby('time.year').sum()
    tmax_occur_weight_mean_yearly = tmax_occur_weight_mean_yearly.squeeze(drop=True)
    return tmax_occur_weight_mean_yearly

def save_result(domain):
    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'

    for f in ['piControl']:
        path = to_file_dir + 'hot_extreme_occur_' + f + '_yearly_all_patterns_' + domain + '.csv'
        if os.path.exists(path) == True:
            pass 
        else:
            data_run_name = []
            for d in dataset_src_run.keys():
                for k in dataset_src_run[d][f]:
                    data_run_name.append(d + '_' + k)
            
            tmax_occur_cmip = pd.DataFrame(columns=data_run_name)
            for d in dataset_src_run.keys():
                for r in dataset_src_run[d][f]:
                    t0 = time_range_tmax[d][f][0]
                    tmax_occur_concat = cal_hot_extreme_occur(dataset_name=d,domain=domain,forcing=f,run=r,time=t0)
                    tmax_occur_concat = tmax_occur_concat.values
                    tmax_occur_concat = pd.Series(tmax_occur_concat)
                    print(tmax_occur_concat)

                    for t in time_range_tmax[d][f]:
                        if t == t0:
                            continue
                        else:
                            tmax_occur = cal_hot_extreme_occur(dataset_name=d,domain=domain,forcing=f,run=r,time=t)
                            tmax_occur = tmax_occur.values
                            tmax_occur = pd.Series(tmax_occur)
                            tmax_occur_concat = pd.concat([tmax_occur_concat,tmax_occur],axis=0,ignore_index=True)
                            print(tmax_occur)
                            print(f + '_' + d + '_' + r + '_' + t)
                    tmax_occur_cmip[d + '_' + r] = tmax_occur_concat
            tmax_occur_cmip.to_csv(to_file_dir + 'hot_extreme_occur_' + f + '_yearly_all_patterns_' + domain + '.csv')
    
save_result('WNA')
save_result('EU')
save_result('EAS')

