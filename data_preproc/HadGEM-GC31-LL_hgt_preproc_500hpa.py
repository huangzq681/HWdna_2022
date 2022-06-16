import xarray as xr
import os

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
sel_level = 500 * 100

var_name = 'zg'
forcs = ['hist-GHG','hist-aer','hist-nat','historical','ssp585']
ensems = ['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3']
starts1 = ['19500101','20000101']
ends1 = ['19991230','20201230']
starts2 = ['19500101','20000101']
ends2 = ['19991230','20141230']
starts3 = ['20500101']
ends3 = ['20991230']

zg0_ds = xr.open_dataset('zg_day_HadGEM3-GC31-LL_ssp585_r3i1p1f3_gn_20150101-20491230.nc')
zg0 = zg0_ds['zg']
mask_lon = (zg0.lon >= lon_min) & (zg0.lon <= lon_max)
mask_lat = (zg0.lat >= lat_min) & (zg0.lat <= lat_max)

def get_file(forc,ensem,start,end):
	file_name = var_name + '_day_HadGEM3-GC31-LL_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '.nc'
	return file_name

def to_file(forc,ensem,start,end):
	file_name = var_name + '_day_HadGEM3-GC31-LL_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '_level' + str(sel_level) + '.nc'
	return file_name

for f in forcs:
	if f in ['hist-GHG','hist-aer','hist-nat']:
		for e in ensems:
			to_file_path = to_file(f,e,starts1[0],ends1[1])
			if os.path.exists(to_file_path):
				continue

			file1 = get_file(f,e,starts1[0],ends1[0])
			zg1_ds = xr.open_dataset(file1)
			zg1 = zg1_ds['zg']
			zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
			zg1_crop_lev = zg1_crop.sel(plev = sel_level)

			file2 = get_file(f,e,starts1[1],ends1[1])
			zg2_ds = xr.open_dataset(file2)
			zg2 = zg2_ds['zg']
			zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
			zg2_crop_lev = zg2_crop.sel(plev = sel_level)

			zg_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
			zg_crop_lev.to_netcdf(to_file_path)
			print(f,e)

	elif f == 'historical':
		for e in ensems:
			to_file_path = to_file(f,e,starts2[0],ends2[1])
			if os.path.exists(to_file_path):
				continue

			file1 = get_file(f,e,starts2[0],ends2[0])
			zg1_ds = xr.open_dataset(file1)
			zg1 = zg1_ds['zg']
			zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
			zg1_crop_lev = zg1_crop.sel(plev = sel_level)

			file2 = get_file(f,e,starts2[1],ends2[1])
			zg2_ds = xr.open_dataset(file2)
			zg2 = zg2_ds['zg']
			zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
			zg2_crop_lev = zg2_crop.sel(plev = sel_level)

			zg_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')

			to_file_path = to_file(f,e,starts2[0],ends2[1])
			zg_crop_lev.to_netcdf(to_file_path)
			print(f,e)

	else:
		for e in ensems:
			to_file_path = to_file(f,e,starts3[0],ends3[0])
			if os.path.exists(to_file_path):
				continue

			file1 = get_file(f,e,starts3[0],ends3[0])
			zg1_ds = xr.open_dataset(file1)
			zg1 = zg1_ds['zg']
			zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
			zg1_crop_lev = zg1_crop.sel(plev = sel_level)

			to_file_path = to_file(f,e,starts3[0],ends3[0])
			zg1_crop_lev.to_netcdf(to_file_path)
			print(f,e)
			






