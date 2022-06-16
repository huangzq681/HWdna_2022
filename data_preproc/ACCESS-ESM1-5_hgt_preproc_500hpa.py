import xarray as xr
import os

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
sel_level = 50000

var_name = 'zg'
forcs = ['ssp585']
ensems = ['r1i1p1f1','r2i1p1f1','r3i1p1f1']
starts3 = ['20500101','20550101','20600101','20650101','20700101','20750101','20800101','20850101','20900101','20950101']
ends3 = ['20541231','20591231','20641231','20691231','20741231','20791231','20841231','20891231','20941231','20991231']

zg0_ds = xr.open_dataset('zg_day_ACCESS-ESM1-5_ssp585_r1i1p1f1_gn_20150101-20191231.nc')
zg0 = zg0_ds['zg']
plev500 = zg0['plev'][3]
mask_lon = (zg0.lon >= lon_min) & (zg0.lon <= lon_max)
mask_lat = (zg0.lat >= lat_min) & (zg0.lat <= lat_max)

def get_file(forc,ensem,start,end):
	file_name = var_name + '_day_ACCESS-ESM1-5_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '.nc'
	return file_name

def to_file(forc,ensem,start,end):
	file_name = var_name + '_day_ACCESS-ESM1-5_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '_level' + str(sel_level) + '.nc'
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
			zg1_crop_lev = zg1_crop.sel(plev = plev500)

			file2 = get_file(f,e,starts1[1],ends1[1])
			zg2_ds = xr.open_dataset(file2)
			zg2 = zg2_ds['zg']
			zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
			zg2_crop_lev = zg2_crop.sel(plev = plev500)

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
			zg1_crop_lev = zg1_crop.sel(plev = plev500)

			file2 = get_file(f,e,starts2[1],ends2[1])
			zg2_ds = xr.open_dataset(file2)
			zg2 = zg2_ds['zg']
			zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
			zg2_crop_lev = zg2_crop.sel(plev = plev500)

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
			zg1_crop_lev = zg1_crop.sel(plev = plev500)

			for i in range(len(starts3)-1):
				file2 = get_file(f,e,starts3[i+1],ends3[i+1])
				zg2_ds = xr.open_dataset(file2)
				zg2 = zg2_ds['zg']
				zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
				zg2_crop_lev = zg2_crop.sel(plev = plev500)

				zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')

			to_file_path = to_file(f,e,starts3[0],ends3[-1])
			zg1_crop_lev.to_netcdf(to_file_path)
			print(f,e)
			






