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
forcs = ['piControl']
ensems = ['r1i1p1f1']
starts1 = [str(i*10) + '0101' for i in range(185,205)]
ends1 = [str(i*10 + 9) + '1231' for i in range(185,205)]

zg0_ds = xr.open_dataset('zg_day_MRI-ESM2-0_piControl_r1i1p1f1_gn_18500101-18591231.nc')
zg0 = zg0_ds['zg']
mask_lon = (zg0.lon >= lon_min) & (zg0.lon <= lon_max)
mask_lat = (zg0.lat >= lat_min) & (zg0.lat <= lat_max)

def get_file(forc,ensem,start,end):
    file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '.nc'
    return file_name

def to_file(forc,ensem,start,end):
    file_name = var_name + '_day_MRI-ESM2-0_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '_level' + str(sel_level) + '.nc'
    return file_name

for f in forcs:
    for e in ensems:
        to_file_path = to_file(f,e,starts1[0],ends1[-1])
        if os.path.exists(to_file_path):
            continue

        file1 = get_file(f,e,starts1[0],ends1[0])
        zg1_ds = xr.open_dataset(file1)
        zg1 = zg1_ds['zg']
        zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
        zg1_crop_lev = zg1_crop.sel(plev = sel_level)

        for t in range(len(starts1)-1):
            file2 = get_file(f,e,starts1[t+1],ends1[t+1])
            zg2_ds = xr.open_dataset(file2)
            zg2 = zg2_ds['zg']
            zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
            zg2_crop_lev = zg2_crop.sel(plev = sel_level)
            zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
        zg1_crop_lev.to_netcdf(to_file_path)
        print(f,e)







