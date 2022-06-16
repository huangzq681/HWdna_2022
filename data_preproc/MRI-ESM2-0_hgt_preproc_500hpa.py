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
ensems = ['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1']
starts1 = ['19500101','19600101','19700101','19800101','19900101','20000101','20100101']
ends1 = ['19591231','19691231','19791231','19891231','19991231','20091231','20191231']
starts2 = ['19500101','19600101','19700101','19800101','19900101','20000101','20100101']
ends2 = ['19591231','19691231','19791231','19891231','19991231','20091231','20141231']
starts3 = ['20550101','20650101','20750101','20850101','20950101']
ends3 = ['20641231','20741231','20841231','20941231','21001231']

zg0_ds = xr.open_dataset('zg_day_MRI-ESM2-0_historical_r1i1p1f1_gn_19500101-19591231.nc')
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
    if f in ['hist-GHG','hist-aer','hist-nat']:
        for e in ensems:
            if e in ['r2i1p1f1','r4i1p1f1']:
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
            else:
                to_file_path = to_file(f,e,starts2[0],ends2[-1])
                if os.path.exists(to_file_path):
                    continue

                file1 = get_file(f,e,starts2[0],ends2[0])
                zg1_ds = xr.open_dataset(file1)
                zg1 = zg1_ds['zg']
                zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
                zg1_crop_lev = zg1_crop.sel(plev = sel_level)

                for t in range(len(starts2)-1):
                    file2 = get_file(f,e,starts2[t+1],ends2[t+1])
                    zg2_ds = xr.open_dataset(file2)
                    zg2 = zg2_ds['zg']
                    zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
                    zg2_crop_lev = zg2_crop.sel(plev = sel_level)
                    zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
                zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)

    elif f == 'historical':
        for e in ensems:
            to_file_path = to_file(f,e,starts2[0],ends2[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts2[0],ends2[0])
            zg1_ds = xr.open_dataset(file1)
            zg1 = zg1_ds['zg']
            zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
            zg1_crop_lev = zg1_crop.sel(plev = sel_level)

            for t in range(len(starts2)-1):
                file2 = get_file(f,e,starts2[t+1],ends2[t+1])
                zg2_ds = xr.open_dataset(file2)
                zg2 = zg2_ds['zg']
                zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
                zg2_crop_lev = zg2_crop.sel(plev = sel_level)
                zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
            zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)

    else:
        for e in ensems:
            to_file_path = to_file(f,e,starts3[0],ends3[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts3[0],ends3[0])
            zg1_ds = xr.open_dataset(file1)
            zg1 = zg1_ds['zg']
            zg1_crop = zg1.where(mask_lon & mask_lat, drop=True)
            zg1_crop_lev = zg1_crop.sel(plev = sel_level)

            for t in range(len(starts3)-1):
                file2 = get_file(f,e,starts3[t+1],ends3[t+1])
                zg2_ds = xr.open_dataset(file2)
                zg2 = zg2_ds['zg']
                zg2_crop = zg2.where(mask_lon & mask_lat, drop=True)
                zg2_crop_lev = zg2_crop.sel(plev = sel_level)
                zg1_crop_lev = xr.concat([zg1_crop_lev,zg2_crop_lev],'time')
            zg1_crop_lev.to_netcdf(to_file_path)
            print(f,e)
            


