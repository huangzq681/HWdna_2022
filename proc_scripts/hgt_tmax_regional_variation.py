import xarray as xr
import numpy as np
import pandas as pd
import os
# import xesmf as xe
# from cartopy.util import add_cyclic_point
from scipy.stats import linregress

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

time_range_hgt = {
    'era5':'1950-2020',
    'jra55':'1958-2014',
    'ncep2':'1979-2020',
    'CanESM5':{
        'historical':'19510101-20141231',
        'hist-GHG':'19510101-20201231',
        'hist-nat':'19510101-20201231',
        'hist-aer':'19510101-20201231',
        'ssp585':'20510101-21001231'
    },
    'HadGEM3-GC31-LL':{
        'historical':'19500101-20141230',
        'hist-GHG':'19500101-20201230',
        'hist-nat':'19500101-20201230',
        'hist-aer':'19500101-20201230',
        'ssp585':'20500101-20991230'
    },
    'MIROC6':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20500101-21001231'
    },
    'IPSL-CM6A-LR':{
        'historical':'19500101-20141231',
        'hist-GHG':'19600101-20201231',
        'hist-nat':'19600101-20201231',
        'hist-aer':'19600101-20201231',
        'ssp585':'20150101-21001231'
    },
    'MRI-ESM2-0':{
        'historical':'19500101-20141231',
        'hist-GHG':'19500101-20141231',
        'hist-nat':'19500101-20141231',
        'hist-aer':'19500101-20141231',
        'ssp585':'20550101-21001231'
    }
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

def sel_domain(dataarray,lon_min,lon_max,lat_min,lat_max):
    mask_lon = (dataarray.lon >= lon_min) & (dataarray.lon <= lon_max)
    mask_lat = (dataarray.lat >= lat_min) & (dataarray.lat <= lat_max)
    dataarray = dataarray.where(mask_lon & mask_lat, drop=True)
    return dataarray

def get_regional_variation(var_name,dataset_name,forcing=None,run=None,domain='NH',method='mean'):
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    var_name_tran = {'hgt':'zg','tmax':'tasmax'}
    if dataset_name in ['era5','jra55','ncep2']:
        filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_hgt[dataset_name] + '.nc'
    else:
        if dataset_name != 'IPSL-CM6A-LR' and var_name == 'hgt':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
        elif dataset_name != 'IPSL-CM6A-LR' and var_name != 'hgt':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_tmax[dataset_name][forcing] + '.nc'
        elif dataset_name == 'IPSL-CM6A-LR' and var_name == 'hgt':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
        elif dataset_name == 'IPSL-CM6A-LR' and var_name != 'hgt':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_tmax[dataset_name][forcing] + '.nc'
        else:
            pass
    if dataset_name in ['era5','jra55','ncep2']:
        data = xr.open_dataarray(filepath)
    else:
        data = xr.open_dataset(filepath)
        data = data[var_name_tran[var_name]]
    if dataset_name == 'era5':
        data = data.rename({'longitude':'lon','latitude':'lat','time':'time'})

    if method == 'mean': 
        data_ds_season_mean = data.resample(time = 'QS-DEC').mean()
    elif method == 'max':
        data_ds_season_mean = data.resample(time = 'QS-DEC').max()
    else:
        pass
    data_ds_summer_mean = data_ds_season_mean.sel(time=data_ds_season_mean['time.season']=='JJA')
    if dataset_name == 'era5' and var_name == 'hgt':
        data_ds_summer_mean = data_ds_summer_mean / 9.80665  ## transform geopotential to geopotential height

    if domain == 'NH':
        pass
    elif domain == 'EAS':
        data_ds_summer_mean = sel_domain(data_ds_summer_mean,lon_min=85,lon_max=125,lat_min=25,lat_max=55)
    elif domain == 'EU':
        data_ds_summer_mean = sel_domain(data_ds_summer_mean,lon_min=10,lon_max=50,lat_min=35,lat_max=65)
    elif domain == 'WNA':
        data_ds_summer_mean = sel_domain(data_ds_summer_mean,lon_min=230,lon_max=258,lat_min=25,lat_max=55)
    else:
        raise Exception('ERROR, domain must be a string object belongs to [\'NH\',\'EAS\',\'EU\',\'WNA\']')

    weights = np.cos(np.deg2rad(data_ds_summer_mean.lat))
    weights.name = "weights"
    data_weighted = data_ds_summer_mean.weighted(weights)
    weighted_mean = data_weighted.mean(("lon", "lat"))
    wm = weighted_mean.values
    wm = pd.DataFrame(wm,index=weighted_mean.time.dt.year,columns=[dataset_name])
    wm.index.name = 'time'
    wm_ano = wm - wm.loc['1979':'2014'].mean()

    return filepath,wm,wm_ano


def save_result(domain='NH'):

    wm_hgt = pd.DataFrame()
    wm_tmax = pd.DataFrame()
    wm_maxTx = pd.DataFrame()
    wm_hgt_ano = pd.DataFrame()
    wm_tmax_ano = pd.DataFrame()
    wm_maxTx_ano = pd.DataFrame()

    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    var_name = 'tmax'
    for d in time_range_tmax.keys():
        if d in ['era5','jra55','ncep2']:
            path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,method='mean',domain=domain)
            if os.path.exists(path) == True:
                print(path)
                wm.columns = [d]
                wm_ano.columns = [d]
                wm_tmax = pd.concat([wm_tmax, wm],axis=1)
                wm_tmax_ano = pd.concat([wm_tmax_ano, wm_ano],axis=1)
            else:
                print('*********'+path)
        else:
            for f in dataset_src_run[d].keys():
                for r in dataset_src_run[d][f]:
                    path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,forcing=f,run=r,method='mean',domain=domain)
                    if os.path.exists(path) == True:
                        print(path)
                        wm.columns = [d + '_' + f + '_' + r]
                        wm_ano.columns = [d + '_' + f + '_' + r]
                        wm_tmax = pd.concat([wm_tmax, wm],axis=1)
                        wm_tmax_ano = pd.concat([wm_tmax_ano, wm_ano],axis=1)
                    else:
                        print('********'+path)
    wm_tmax.to_csv(to_file_dir + var_name + 'mean_regional_variation_' + domain + '.csv')
    wm_tmax_ano.to_csv(to_file_dir + var_name + 'mean_anomaly_regional_variation_' + domain + '.csv')

    var_name = 'tmax'
    for d in time_range_tmax.keys():
        if d in ['era5','jra55','ncep2']:
            path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,method='max',domain=domain)
            if os.path.exists(path) == True:
                print(path)
                wm.columns = [d]
                wm_ano.columns = [d]
                wm_maxTx = pd.concat([wm_maxTx, wm],axis=1)
                wm_maxTx_ano = pd.concat([wm_maxTx_ano, wm_ano],axis=1)
            else:
                print('*********'+path)
        else:
            for f in dataset_src_run[d].keys():
                for r in dataset_src_run[d][f]:
                    path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,forcing=f,run=r,method='max',domain=domain)
                    if os.path.exists(path) == True:
                        print(path)
                        wm.columns = [d + '_' + f + '_' + r]
                        wm_ano.columns = [d + '_' + f + '_' + r]
                        wm_maxTx = pd.concat([wm_maxTx, wm],axis=1)
                        wm_maxTx_ano = pd.concat([wm_maxTx_ano, wm_ano],axis=1)
                    else:
                        print('********'+path)
    wm_maxTx.to_csv(to_file_dir + var_name + 'max_regional_variation_' + domain + '.csv')
    wm_maxTx_ano.to_csv(to_file_dir + var_name + 'max_anomaly_regional_variation_' + domain + '.csv')

    var_name = 'hgt'
    for d in time_range_hgt.keys():
        if d in ['era5','jra55','ncep2']:
            path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,method='mean',domain=domain)
            if os.path.exists(path) == True:
                print(path)
                wm.columns = [d]
                wm_ano.columns = [d]
                wm_hgt = pd.concat([wm_hgt, wm],axis=1)
                wm_hgt_ano = pd.concat([wm_hgt_ano, wm_ano],axis=1)
            else:
                print('*********'+path)
        else:
            for f in dataset_src_run[d].keys():
                for r in dataset_src_run[d][f]:
                    path, wm, wm_ano = get_regional_variation(var_name,dataset_name=d,forcing=f,run=r,method='mean',domain=domain)
                    if os.path.exists(path) == True:
                        print(path)
                        wm.columns = [d + '_' + f + '_' + r]
                        wm_ano.columns = [d + '_' + f + '_' + r]
                        wm_hgt = pd.concat([wm_hgt, wm],axis=1)
                        wm_hgt_ano = pd.concat([wm_hgt_ano, wm_ano],axis=1)
                    else:
                        print('********'+path)
    wm_hgt.to_csv(to_file_dir + var_name + 'mean_regional_variation_' + domain + '.csv')
    wm_hgt_ano.to_csv(to_file_dir + var_name + 'mean_anomaly_regional_variation_' + domain + '.csv')

save_result(domain='NH')
save_result(domain='EAS')
save_result(domain='WNA')
save_result(domain='EU')