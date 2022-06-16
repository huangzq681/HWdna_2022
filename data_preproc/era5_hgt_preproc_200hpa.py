import xarray as xr

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
sel_level = 200
year_start = 1950
year_end = 2020

i0 = year_start
hgt_i_ds = xr.open_dataset('geopotential_daily_era5_'+str(i0)+'.nc')
hgt_i = hgt_i_ds['z']
mask_lon = (hgt_i.longitude >= lon_min) & (hgt_i.longitude <= lon_max)
mask_lat = (hgt_i.latitude >= lat_min) & (hgt_i.latitude <= lat_max)
hgt_i_crop = hgt_i.where(mask_lon & mask_lat, drop=True)
hgt_all_crop_500lev = hgt_i_crop.sel(level=sel_level)

for i in range(year_start+1,year_end+1):
	hgt_i_ds = xr.open_dataset('geopotential_daily_era5_'+str(i)+'.nc')
	hgt_i = hgt_i_ds['z']
	hgt_i_crop = hgt_i.where(mask_lon & mask_lat, drop=True)
	hgt_i_crop_500lev = hgt_i_crop.sel(level=sel_level)

	hgt_all_crop_500lev = xr.concat([hgt_all_crop_500lev,hgt_i_crop_500lev],'time')
	print(i)

to_path = 'geopotential_x' + str(lon_min) + '-' + str(lon_max) + '_y' + str(lat_min) + '-' + str(lat_max) + '_t' + str(year_start) + '-' + str(year_end) + '_' + 'level' + str(sel_level) + '.nc'
hgt_all_crop_500lev.to_netcdf(to_path)
