import xarray as xr

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
year_start = 1950
year_end = 2020

i0 = year_start
mx2t_i_ds = xr.open_dataset('maximum_2m_temperature_since_previous_post_processing_daily_era5_'+str(i0)+'.nc')
mx2t_i = mx2t_i_ds['mx2t']
mask_lon = (mx2t_i.longitude >= lon_min) & (mx2t_i.longitude <= lon_max)
mask_lat = (mx2t_i.latitude >= lat_min) & (mx2t_i.latitude <= lat_max)
mx2t_i_crop_0 = mx2t_i.where(mask_lon & mask_lat, drop=True)

for i in range(year_start+1,year_end+1):
	mx2t_i_ds = xr.open_dataset('maximum_2m_temperature_since_previous_post_processing_daily_era5_'+str(i)+'.nc')
	mx2t_i = mx2t_i_ds['mx2t']
	mx2t_i_crop = mx2t_i.where(mask_lon & mask_lat, drop=True)

	mx2t_i_crop_0 = xr.concat([mx2t_i_crop_0,mx2t_i_crop],'time')
	print(i)

to_path = 'mx2t_x' + str(lon_min) + '-' + str(lon_max) + '_y' + str(lat_min) + '-' + str(lat_max) + '_t' + str(year_start) + '-' + str(year_end) + '.nc'
mx2t_i_crop_0.to_netcdf(to_path)
