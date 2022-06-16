'''
This script is for determining the winner of JJA 500hpa geopotential height anomalies relative to multi-reanalysis mean for era5, ncep2, jra55 and
CMIP6 historical dataset by computing the euclidean distance. 
And then the occurrence of each SOM node are computed for the dataset.  
'''
import xarray as xr
from minisom import MiniSom
import pandas as pd
import numpy as np
from scipy.spatial import distance
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

dataset_src_run = {
    'CanESM5':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1']
    },
    'HadGEM3-GC31-LL':{
        'historical':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3']
    },
    'MIROC6':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1']
    },
    'IPSL-CM6A-LR':{
        'historical':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1']
    },
    'MRI-ESM2-0':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']
    },
}

time_range_hgt = {
    'era5':'1950-2020',
    'jra55':'1958-2014',
    'ncep2':'1979-2020',
    'CanESM5':{
        'historical':'19510101-20141231'
    },
    'HadGEM3-GC31-LL':{
        'historical':'19500101-20141230'
    },
    'MIROC6':{
        'historical':'19500101-20141231'
    },
    'IPSL-CM6A-LR':{
        'historical':'19500101-20141231'
    },
    'MRI-ESM2-0':{
        'historical':'19500101-20141231'
    }
}

domain_lonlat = {
    'EAS':{'lon_min':90,'lon_max':130,'lat_min':30,'lat_max':60},
    'EU':{'lon_min':10,'lon_max':50,'lat_min':35,'lat_max':65},
    # 'WNA':{'lon_min':230,'lon_max':262,'lat_min':30,'lat_max':54},
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
    var_name_tran = {'hgt':'zg'}
    if dataset_name in ['era5','jra55','ncep2']:
        filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_hgt[dataset_name] + '.nc'
    else:
        if dataset_name != 'IPSL-CM6A-LR':
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gn_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
        else:
            filepath = filedir + dataset_name + '/' + var_name_tran[var_name] + '_day_' + dataset_name + '_' + forcing + '_' + run + '_gr_' + time_range_hgt[dataset_name][forcing] + '_level50000.nc'
        
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
    if dataset_name == 'era5' and var_name == 'hgt':
        data_summer = data_summer / 9.80665  ## transform geopotential to geopotential height
    # data_summer = data_summer.interp(lat=target_griddes['lat'],lon=target_griddes['lon'],method='nearest')

    if domain == 'EAS':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['EAS']['lon_min'],lon_max=domain_lonlat['EAS']['lon_max'],lat_min=domain_lonlat['EAS']['lat_min'],lat_max=domain_lonlat['EAS']['lat_max'])
    elif domain == 'EU':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['EU']['lon_min'],lon_max=domain_lonlat['EU']['lon_max'],lat_min=domain_lonlat['EU']['lat_min'],lat_max=domain_lonlat['EU']['lat_max'])
    elif domain == 'WNA':
        data_summer = sel_lonlat_range(data_summer,lon_min=domain_lonlat['WNA']['lon_min'],lon_max=domain_lonlat['WNA']['lon_max'],lat_min=domain_lonlat['WNA']['lat_min'],lat_max=domain_lonlat['WNA']['lat_max'])
    else:
        raise Exception('ERROR, domain must be a string object belongs to [\'EAS\',\'EU\',\'WNA\']')

    return data_summer

def determine_winner(regrid_domain_dat,domain):

    hgt_ds_summer = regrid_domain_dat
    hgt_ds_summer = hgt_ds_summer.groupby('time.dayofyear') - hgt_ds_summer.groupby('time.dayofyear').mean('time')

    # regrid
    hgt_ds_summer = hgt_ds_summer.interp(lat=target_griddes_sub[domain]['lat'],lon=target_griddes_sub[domain]['lon'])
    hgt_ds_summer = hgt_ds_summer.interpolate_na(dim='lon',method='linear',fill_value='extrapolate')
    hgt_ds_summer = hgt_ds_summer.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')

    hgt_ds_summer_stack = hgt_ds_summer.stack(z=('lon','lat'))
    hgt_ds_summer_stack = hgt_ds_summer_stack.values
    hgt_ds_summer_stack = pd.DataFrame(hgt_ds_summer_stack)
    hgt_ds_summer_stack = hgt_ds_summer_stack.to_numpy()
 
    som_weight_domain_path = '/lustre07/scratch/xtan/hzq/HWdna/procData/som_pattern_reanalyses_avg_weights_' + domain + '.csv'
    som_weight_domain = pd.read_csv(som_weight_domain_path,index_col=0)
    winners = []

    for i in range(len(hgt_ds_summer_stack)):
        dist1 = distance.euclidean(hgt_ds_summer_stack[i],som_weight_domain['grid:0-0'])
        dist2 = distance.euclidean(hgt_ds_summer_stack[i],som_weight_domain['grid:1-0'])
        dist3 = distance.euclidean(hgt_ds_summer_stack[i],som_weight_domain['grid:2-0'])
        dist4 = distance.euclidean(hgt_ds_summer_stack[i],som_weight_domain['grid:3-0'])
        dists = np.array([dist1,dist2,dist3,dist4])
        winner = dists.argmin()
        winners.append(winner)

    return winners

def save_result(domain):

    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    som_winner_df = pd.DataFrame()
    time_range = pd.date_range('1979-01-01','2014-12-31')
    time_range = time_range[time_range.month.isin([6,7,8])]
    time_range2 = time_range[time_range.day != 31]

    for d in time_range_hgt.keys():
        if d in ['era5','jra55','ncep2']:
            rd_dat = sel_domain(d,domain)
            som_winner = determine_winner(rd_dat,domain)
            som_winner = pd.Series(som_winner,index=time_range)
            som_winner_df[d] = som_winner
            print(d)
        else:
            for f in dataset_src_run[d].keys():
                for r in dataset_src_run[d][f]:
                    rd_dat = sel_domain(dataset_name=d,domain=domain,forcing=f,run=r)
                    som_winner = determine_winner(rd_dat,domain)
                    if d != 'HadGEM3-GC31-LL':
                        som_winner = pd.Series(som_winner,index=time_range)
                    else:
                        som_winner = pd.Series(som_winner,index=time_range2)
                    som_winner_df[d+'_'+r] = som_winner
                    print(d + '_' + f + '_' + r)
    
    som_winner_df.to_csv(to_file_dir + 'som_winner_relative_to_reanalyses_mean_historical_' + domain + '.csv')

save_result('EAS')
save_result('EU')
save_result('WNA')
