'''
This script is for calculating 4-grid SOM for 3 regions of interest in the CMIP6 ALL forcing dataset, 
and determine the similarities between models and reanalysis.
Fisrtly, we apply the SOM analysis for the 3 regions using the reanalysis and CMIP6 dataset based on their native spatial resolution.
Then, we regrid all the output patterns and calculate the pattern correlations between reanalyses and CMIP6 models for the historical period 1979~2014. 
Since the outputs of SOM analysis are random arranged, we determine the pattern correlations between dataset according to the SOM pairs. 
Finally, a boxplot is drawn to illustrate the performance of CMIP6 models to simulate the historical circulation anomalies.
'''
import xarray as xr
from minisom import MiniSom
import pandas as pd
import numpy as np
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
        raise Exception('ERROR, domain must be a string object belongs to [\'NH\',\'EAS\',\'EU\',\'WNA\']')

    return data_summer

## do m_neurons x n_neurons SOM analysis
'''
The 500hpa geopotential height anomalies in the Observed and ALL forcing data are computed by subtracting the seasonal cycle (calendar-day mean) 
of the historical period (1979~2014) from each grid cell as D.E Horton et al do. The anomalies from different external forcings are computed relatived to 
the historical period (1979~2014) of ALL forcing at the same way to generate comparative data.

ref: Horton, Daniel E., Nathaniel C. Johnson, Deepti Singh, Daniel L. Swain, Bala Rajaratnam, and Noah S. Diffenbaugh. 
“Contribution of Changes in Atmospheric Circulation Patterns to Extreme Temperature Trends.” Nature 522, no. 7557 (June 2015): 465–69. https://doi.org/10/f7hcw7.
'''
def som_analysis(regrid_domain_dat,n_neurons,m_neurons,domain):

    hgt_ds_summer = regrid_domain_dat
    lon0 = hgt_ds_summer.coords['lon']
    lat0 = hgt_ds_summer.coords['lat']
    hgt_ds_summer = hgt_ds_summer.groupby('time.dayofyear') - hgt_ds_summer.groupby('time.dayofyear').mean('time')

    hgt_ds_summer_stack = hgt_ds_summer.stack(z=('lon','lat'))
    hgt_ds_summer_stack = hgt_ds_summer_stack.values
    hgt_ds_summer_stack = pd.DataFrame(hgt_ds_summer_stack)
    hgt_ds_summer_stack = hgt_ds_summer_stack.to_numpy()
 
    som = MiniSom(n_neurons, m_neurons, hgt_ds_summer_stack.shape[1], sigma=0.5, learning_rate=.5,
                neighborhood_function='gaussian', random_seed=0)
    # som.pca_weights_init(hgt_ds_summer_stack)
    som.train(hgt_ds_summer_stack, 10000) # trains the SOM with 10000 iterations

    som_grid_df = pd.DataFrame()
    for n in range(n_neurons):
        for m in range(m_neurons):
            col_name = 'grid:' + str(n) + '-' + str(m)
            som_grid = som.get_weights()[n,m,:]
            som_grid = som_grid.reshape([len(lat0),len(lon0)],order='F')
            som_grid = xr.DataArray(
                data=som_grid,dims=['lat','lon'],
                coords=dict(
                    lon = lon0, lat = lat0
                ))
            som_grid = som_grid.interp(lat=target_griddes_sub[domain]['lat'],lon=target_griddes_sub[domain]['lon'])
            som_grid = som_grid.interpolate_na(dim='lon',method='linear',fill_value='extrapolate')
            som_grid = som_grid.interpolate_na(dim='lat',method='linear',fill_value='extrapolate')
            som_grid_stack = som_grid.stack(z=('lon','lat'))
            som_grid_df[col_name] = som_grid_stack      
    
    winner_coordinates = np.array([som.winner(x) for x in hgt_ds_summer_stack]).T
    cluster_index = np.ravel_multi_index(winner_coordinates, (n_neurons,m_neurons))

    return cluster_index, som_grid_df

def save_result(domain,n_neurons,m_neurons):

    to_file_dir = '/home/xtan/scratch/hzq/HWdna/procData/'
    som_pattern_grid_writer = pd.ExcelWriter(to_file_dir + 'som_pattern_grid_' + str(n_neurons) + '-' + str(m_neurons) + '_historical_'+domain+'.xlsx', engine='xlsxwriter')
    som_winner_df = pd.DataFrame()
    time_range = pd.date_range('1979-01-01','2014-12-31')
    time_range = time_range[time_range.month.isin([6,7,8])]
    time_range2 = time_range[time_range.day != 31]

    for d in time_range_hgt.keys():
        if d in ['era5','jra55','ncep2']:
            rd_dat = sel_domain(d,domain)
            som_winner,som_grid = som_analysis(rd_dat,n_neurons,m_neurons,domain)
            som_winner = pd.Series(som_winner,index=time_range)
            som_winner_df[d] = som_winner
            som_grid.to_excel(som_pattern_grid_writer, sheet_name=d)
            print(d)
        else:
            for f in dataset_src_run[d].keys():
                for r in dataset_src_run[d][f]:
                    rd_dat = sel_domain(dataset_name=d,domain=domain,forcing=f,run=r)
                    som_winner,som_grid = som_analysis(rd_dat,n_neurons,m_neurons,domain)
                    if d != 'HadGEM3-GC31-LL':
                        som_winner = pd.Series(som_winner,index=time_range)
                    else:
                        som_winner = pd.Series(som_winner,index=time_range2)
                    som_winner_df[d+'_'+r] = som_winner
                    som_grid.to_excel(som_pattern_grid_writer, sheet_name=d+'_'+r)
                    print(d + '_' + f + '_' + r)
    
    som_winner_df.to_csv(to_file_dir + 'som_winner_historical_' + domain + '.csv')
    som_pattern_grid_writer.save()

# save_result('EAS',4,1)
# save_result('EU',4,1)
save_result('WNA',4,1)
