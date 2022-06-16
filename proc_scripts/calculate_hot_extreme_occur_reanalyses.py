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
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)

time_range_tmax = {
    'jra55':'1958-2014',
    'era5':'1950-2020',
    'ncep2':'1979-2020'}

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

def sel_domain(dataset_name,domain):
    var_name = 'tmax'
    filedir = '/home/xtan/scratch/hzq/HWdna/rawData/'
    filepath = filedir + dataset_name + '/'+ dataset_name +'_daily_' + var_name + '_' + time_range_tmax[dataset_name] + '.nc'
    data = xr.open_dataarray(filepath)
    if 'height' in data.dims:
        data = data.squeeze(drop=True)
    if dataset_name == 'era5':
        data = data.rename({'longitude':'lon','latitude':'lat','time':'time'})
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


def cal_hot_extreme_occur(dataset_name,domain,forcing=None,run=None):
    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    rd_dat = sel_domain(dataset_name=dataset_name,domain=domain)

    season_cycle_quantile_path = to_file_dir + dataset_name + '_tmax_season_cycle_quantile_1979-2014' + domain + '.nc'
    season_cycle_quantile = xr.open_dataarray(season_cycle_quantile_path)
    season_cycle_quantile.squeeze(drop=True)

    tmax_ano = cal_hot_extreme_anomalies(regrid_domain_dat = rd_dat, season_cycle_data = season_cycle_quantile, domain = domain)
    tmax_ano = tmax_ano.squeeze(drop=True)

    som_winner = pd.read_csv('/home/xtan/scratch/hzq/HWdna/procData/' + 'som_winner_relative_to_reanalyses_mean_historical_' + domain + '.csv', index_col=0)

    years = [i[:4] for i in som_winner.index]
    som_winner_1 = som_winner[dataset_name] == 0
    som_winner_2 = som_winner[dataset_name] == 1
    som_winner_3 = som_winner[dataset_name] == 2
    som_winner_4 = som_winner[dataset_name] == 3

    som_winner_1_np = som_winner_1.to_numpy() * 1
    som_winner_2_np = som_winner_2.to_numpy() * 1
    som_winner_3_np = som_winner_3.to_numpy() * 1
    som_winner_4_np = som_winner_4.to_numpy() * 1

    tmax_ano_patt_1 = tmax_ano * som_winner_1_np[:,np.newaxis,np.newaxis]
    tmax_ano_patt_2 = tmax_ano * som_winner_2_np[:,np.newaxis,np.newaxis]
    tmax_ano_patt_3 = tmax_ano * som_winner_3_np[:,np.newaxis,np.newaxis]
    tmax_ano_patt_4 = tmax_ano * som_winner_4_np[:,np.newaxis,np.newaxis]

    tmax_ano_patt_ishot_1 = (tmax_ano_patt_1 > 0) * 1
    tmax_ano_patt_ishot_2 = (tmax_ano_patt_2 > 0) * 1
    tmax_ano_patt_ishot_3 = (tmax_ano_patt_3 > 0) * 1
    tmax_ano_patt_ishot_4 = (tmax_ano_patt_4 > 0) * 1

    tmax_ano_patt_1_occur = tmax_ano_patt_ishot_1.groupby('time.year').sum()
    tmax_ano_patt_2_occur = tmax_ano_patt_ishot_2.groupby('time.year').sum()
    tmax_ano_patt_3_occur = tmax_ano_patt_ishot_3.groupby('time.year').sum()
    tmax_ano_patt_4_occur = tmax_ano_patt_ishot_4.groupby('time.year').sum()
    
    tmax_ano_patt_1_occur.to_netcdf(to_file_dir + dataset_name + '_' + 'hot_extreme_patt1_occur_' + domain + '.nc')
    tmax_ano_patt_2_occur.to_netcdf(to_file_dir + dataset_name + '_' + 'hot_extreme_patt2_occur_' + domain + '.nc')
    tmax_ano_patt_3_occur.to_netcdf(to_file_dir + dataset_name + '_' + 'hot_extreme_patt3_occur_' + domain + '.nc')
    tmax_ano_patt_4_occur.to_netcdf(to_file_dir + dataset_name + '_' + 'hot_extreme_patt4_occur_' + domain + '.nc')
    
    weights = np.cos(np.deg2rad(tmax_ano_patt_1.lat))
    weights.name = "weights"
    
    tmax_area_weight_int_1 = tmax_ano_patt_1.weighted(weights)
    tmax_area_weight_int_2 = tmax_ano_patt_2.weighted(weights)
    tmax_area_weight_int_3 = tmax_ano_patt_3.weighted(weights)
    tmax_area_weight_int_4 = tmax_ano_patt_4.weighted(weights)

    tmax_area_weight_mean_int_1 = tmax_area_weight_int_1.mean(("lon", "lat"))
    tmax_area_weight_mean_int_2 = tmax_area_weight_int_2.mean(("lon", "lat"))
    tmax_area_weight_mean_int_3 = tmax_area_weight_int_3.mean(("lon", "lat"))
    tmax_area_weight_mean_int_4 = tmax_area_weight_int_4.mean(("lon", "lat"))
    
    tmax_area_weight_mean_int_yearly_1 = tmax_area_weight_mean_int_1.groupby('time.year').sum()
    tmax_area_weight_mean_int_yearly_2 = tmax_area_weight_mean_int_2.groupby('time.year').sum()
    tmax_area_weight_mean_int_yearly_3 = tmax_area_weight_mean_int_3.groupby('time.year').sum()
    tmax_area_weight_mean_int_yearly_4 = tmax_area_weight_mean_int_4.groupby('time.year').sum()

    patt_occur_yearly_1 = som_winner_1.groupby(years).sum()
    patt_occur_yearly_2 = som_winner_2.groupby(years).sum()
    patt_occur_yearly_3 = som_winner_3.groupby(years).sum()
    patt_occur_yearly_4 = som_winner_4.groupby(years).sum()
    
    tmax_area_weight_mean_int_yearly_1 = tmax_area_weight_mean_int_yearly_1 / patt_occur_yearly_1
    tmax_area_weight_mean_int_yearly_2 = tmax_area_weight_mean_int_yearly_2 / patt_occur_yearly_2
    tmax_area_weight_mean_int_yearly_3 = tmax_area_weight_mean_int_yearly_3 / patt_occur_yearly_3
    tmax_area_weight_mean_int_yearly_4 = tmax_area_weight_mean_int_yearly_4 / patt_occur_yearly_4
    
    tmax_area_weight_mean_int_yearly = pd.DataFrame(
        {'Pattern1':tmax_area_weight_mean_int_yearly_1,'Pattern2':tmax_area_weight_mean_int_yearly_2,'Pattern3':tmax_area_weight_mean_int_yearly_3,'Pattern4':tmax_area_weight_mean_int_yearly_4},columns=['Pattern1','Pattern2','Pattern3','Pattern4'])
    tmax_area_weight_mean_int_yearly.to_csv(to_file_dir + dataset_name + '_' + 'hot_extreme_per_pattern_occur_' + domain + '.csv')

def save_result(domain):
    for d in time_range_tmax.keys():
        if d in ['era5','jra55','ncep2']:
            cal_hot_extreme_occur(dataset_name=d,domain=domain)
            print(d)

# save_result('WNA')
# save_result('EU')
# save_result('EAS')

for d in ['era5','jra55','ncep2']:
    for region in ['EU','EAS','WNA']:
        cal_hot_extreme_occur(dataset_name=d,domain=region)
        cal_hot_extreme_occur(dataset_name=d,domain=region)
        cal_hot_extreme_occur(dataset_name=d,domain=region)
        print(d + '_' + region)

