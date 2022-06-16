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
ensems1 = ['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1']
ensems2 = ['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r6i1p1f1','r14i1p1f1']

starts1 = ['19600101','20100101']
ends1 = ['20091231','20201231']
starts2 = ['19500101']
ends2 = ['20141231']
starts3 = ['20150101']
ends3 = ['21001231']

zg0_ds = xr.open_dataset('zg_day_IPSL-CM6A-LR_historical_r1i1p1f1_gr_19500101-20141231.nc')
zg0 = zg0_ds['zg']
mask_lon = (zg0.lon >= lon_min) & (zg0.lon <= lon_max)
mask_lat = (zg0.lat >= lat_min) & (zg0.lat <= lat_max)

def get_file(forc,ensem,start,end):
    if forc in ['hist-GHG','hist-aer','hist-nat']:
        file_name = var_name + '_Eday_IPSL-CM6A-LR_' + forc + '_' + ensem + '_gr_' + start + '-' + end + '.nc'    
    else:
        file_name = var_name + '_day_IPSL-CM6A-LR_' + forc + '_' + ensem + '_gr_' + start + '-' + end + '.nc'
    return file_name

def to_file(forc,ensem,start,end):
    file_name = var_name + '_day_IPSL-CM6A-LR_' + forc + '_' + ensem + '_gr_' + start + '-' + end + '_level' + str(sel_level) + '.nc'
    return file_name

for f in forcs:

    if f in ['hist-GHG','hist-aer','hist-nat']:
        for e in ensems1:
            to_file_path = to_file(f,e,starts1[0],ends1[-1])
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
            zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
            zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)

    elif f == 'historical':
        for e in ensems1:
            to_file_path = to_file(f,e,starts2[0],ends2[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts2[0],ends2[0])
            zg1_ds = xr.open_dataset(file1)
            zg1 = zg1_ds['zg']
            zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
            zg1_crop_lev = zg1_crop.sel(plev = sel_level)
            zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)

    else:
        for e in ensems2:
            to_file_path = to_file(f,e,starts3[0],ends3[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts3[0],ends3[0])
            zg1_ds = xr.open_dataset(file1)
            zg1 = zg1_ds['zg']
            zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
            zg1_crop_lev = zg1_crop.sel(plev = sel_level)
            zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)
            






