import xarray as xr

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
year_start = 1979
year_end = 2020

i0 = year_start
tmax_i_ds = xr.open_dataset('tmax.2m.gauss.'+str(i0)+'.nc')
tmax_i = tmax_i_ds['tmax']
mask_lon = (tmax_i.lon >= lon_min) & (tmax_i.lon <= lon_max)
mask_lat = (tmax_i.lat >= lat_min) & (tmax_i.lat <= lat_max)
tmax_i_crop_0 = tmax_i.where(mask_lon & mask_lat, drop=True)

for i in range(year_start+1,year_end+1):
	tmax_i_ds = xr.open_dataset('tmax.2m.gauss.'+str(i)+'.nc')
	tmax_i = tmax_i_ds['tmax']
	tmax_i_crop = tmax_i.where(mask_lon & mask_lat, drop=True)
	tmax_i_crop_0 = xr.concat([tmax_i_crop_0,tmax_i_crop],'time')
	print(i)

to_path = 'tmax_x' + str(lon_min) + '-' + str(lon_max) + '_y' + str(lat_min) + '-' + str(lat_max) + '_t' + str(year_start) + '-' + str(year_end)  + '.nc'
tmax_i_crop_0.to_netcdf(to_path)
