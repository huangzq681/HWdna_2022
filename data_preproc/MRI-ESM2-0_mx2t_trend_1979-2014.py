import xarray as xr
from scipy.stats import linregress

time_start = '1979-01-01'
time_end   = '2014-12-31'

## select level and date range
var_name = 'tasmax'
forcs = ['hist-GHG','hist-aer','hist-nat','historical']
ensems = ['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']

def _compute_slope(var):
    """
    Private function to compute slopes at each grid cell using
    linregress. 
    """
    # x = np.arange(len(y))
    # slp = np.polyfit(x, y, 1)[0]
    slp = linregress(range(len(var)),var).slope
    return slp

def _compute_sig(var):
    """
    Private function to compute slopes at each grid cell using
    linregress. 
    """
    # slp = linregress(range(len(var)),var).slope
    sig = linregress(range(len(var)),var).pvalue
    return sig

def trend_cal(data):
    """
    Computes linear slope (m) at each grid cell.
    
    Args:
        da: xarray DataArray to compute slopes for
        
    Returns:
        xarray DataArray with slopes computed at each grid cell.
    """
    # apply_ufunc can apply a raw numpy function to a grid.
    # 
    # vectorize is only needed for functions that aren't already
    # vectorized. You don't need it for polyfit in theory, but it's
    # good to use when using things like np.cov.
    #
    # dask='parallelized' parallelizes this across dask chunks. It requires
    # an output_dtypes of the numpy array datatype coming out.
    #
    # input_core_dims should pass the dimension that is being *reduced* by this operation,
    # if one is being reduced.
    slopes = xr.apply_ufunc(_compute_slope,
                            data,
                            vectorize=True,
                            dask='parallelized', 
                            input_core_dims=[['time']],
                            output_dtypes=[float],
                            )
    sig = xr.apply_ufunc(_compute_sig,
                            data,
                            vectorize=True,
                            dask='parallelized', 
                            input_core_dims=[['time']],
                            output_dtypes=[float],
                            )
    return slopes,sig

def get_file(forc,ensem,start='19500101',end='20201231'):
    if forc != 'historical':
        if ensem in ['r2i1p1f1','r4i1p1f1']:
            file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + end +'.nc'
        else:
            file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + '20141231' +'.nc'
    else:
        file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + '20141231' + '.nc'
    return file_name

def to_file(forc,ensem,content,start=time_start,end=time_end):
    file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '_' + content + '.nc'
    return file_name

for f in forcs:
    for e in ensems:
        mx2t_path = get_file(f,e)
        mx2t_ds = xr.open_dataset(mx2t_path)
        mx2t_ds = mx2t_ds.sel(time=slice(time_start,time_end))
        mx2t_ds_season_mean = mx2t_ds.resample(time = 'QS-DEC').mean()
        mx2t_ds_summer_mean = mx2t_ds_season_mean.sel(time=mx2t_ds_season_mean['time.season']=='JJA')
        maxTX_ds_season = mx2t_ds.resample(time = 'QS-DEC').max() ## maximum of daily maximum temperature
        maxTX_ds_summer = maxTX_ds_season.sel(time=maxTX_ds_season['time.season']=='JJA')

        mx2t_ds_summer_mean_clim = mx2t_ds_summer_mean.mean(dim='time')
        mx2t_ds_summer_mean_trend, mx2t_ds_summer_mean_sig = trend_cal(mx2t_ds_summer_mean)
        mx2t_ds_summer_mean_clim.to_netcdf(to_file(f,e,'seasonmean_clim'))
        mx2t_ds_summer_mean_trend.to_netcdf(to_file(f,e,'seasonmean_trend'))
        mx2t_ds_summer_mean_sig.to_netcdf(to_file(f,e,'seasonmean_sig'))

        maxTX_ds_summer_clim = maxTX_ds_summer.mean(dim='time')
        maxTX_ds_summer = maxTX_ds_summer['tasmax']
        maxTX_ds_summer_trend, maxTX_ds_summer_sig = trend_cal(maxTX_ds_summer)
        maxTX_ds_summer_clim.to_netcdf(to_file(f,e,'seasonmax_clim'))
        maxTX_ds_summer_trend.to_netcdf(to_file(f,e,'seasonmax_trend'))
        maxTX_ds_summer_sig.to_netcdf(to_file(f,e,'seasonmax_sig'))

        print(f,e)

        
