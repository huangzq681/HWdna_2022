import xarray as xr
import os

## select roi range
lon_min = 0
lon_max = 360
lat_min = 0
lat_max = 90

## select level and date range
var_name = 'tasmax'
forcs = ['hist-GHG','hist-aer','hist-nat','historical','ssp585']
ensems = ['r1i1p1f1','r2i1p1f1','r3i1p1f1']
starts1 = ['19500101','19600101','19700101','19800101','19900101','20000101','20100101']
ends1 = ['19591231','19691231','19791231','19891231','19991231','20091231','20191231']
starts2 = ['19500101','19600101','19700101','19800101','19900101','20000101','20100101']
ends2 = ['19591231','19691231','19791231','19891231','19991231','20091231','20141231']
starts3 = ['20550101','20650101','20750101','20850101','20950101']
ends3 = ['20641231','20741231','20841231','20941231','21001231']

tasmax0_ds = xr.open_dataset('tasmax_day_MIROC6_ssp585_r1i1p1f1_gn_20150101-20241231.nc')
tasmax0 = tasmax0_ds['tasmax']
mask_lon = (tasmax0.lon >= lon_min) & (tasmax0.lon <= lon_max)
mask_lat = (tasmax0.lat >= lat_min) & (tasmax0.lat <= lat_max)

def get_file(forc,ensem,start,end):
    file_name = var_name + '_day_MIROC6_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '.nc'
    return file_name

def to_file(forc,ensem,start,end):
    file_name = var_name + '_day_MIROC6_' + forc + '_' + ensem + '_gn_' + start + '-' + end + '.nc'
    return file_name

for f in forcs:
    if f in ['hist-GHG','hist-aer','hist-nat']:
        for e in ensems:
            to_file_path = to_file(f,e,starts1[0],ends1[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts1[0],ends1[0])
            tasmax1_ds = xr.open_dataset(file1)
            tasmax1 = tasmax1_ds['tasmax']
            tasmax1_crop = tasmax1.where(mask_lon & mask_lat, drop=True)

            for t in range(len(starts1)-1):
                file2 = get_file(f,e,starts1[t+1],ends1[t+1])
                tasmax2_ds = xr.open_dataset(file2)
                tasmax2 = tasmax2_ds['tasmax']
                tasmax2_crop = tasmax2.where(mask_lon & mask_lat, drop=True)
                tasmax1_crop = xr.concat([tasmax1_crop,tasmax2_crop],'time')
            tasmax1_crop.to_netcdf(to_file_path)
            print(f,e)

    elif f == 'historical':
        for e in ensems:
            to_file_path = to_file(f,e,starts2[0],ends2[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts2[0],ends2[0])
            tasmax1_ds = xr.open_dataset(file1)
            tasmax1 = tasmax1_ds['tasmax']
            tasmax1_crop = tasmax1.where(mask_lon & mask_lat, drop=True)

            for t in range(len(starts2)-1):
                file2 = get_file(f,e,starts2[t+1],ends2[t+1])
                tasmax2_ds = xr.open_dataset(file2)
                tasmax2 = tasmax2_ds['tasmax']
                tasmax2_crop = tasmax2.where(mask_lon & mask_lat, drop=True)
                tasmax1_crop = xr.concat([tasmax1_crop,tasmax2_crop],'time')
            tasmax1_crop.to_netcdf(to_file_path)
            print(f,e)

    else:
        for e in ensems:
            to_file_path = to_file(f,e,starts3[0],ends3[-1])
            if os.path.exists(to_file_path):
                continue

            file1 = get_file(f,e,starts3[0],ends3[0])
            tasmax1_ds = xr.open_dataset(file1)
            tasmax1 = tasmax1_ds['tasmax']
            tasmax1_crop = tasmax1.where(mask_lon & mask_lat, drop=True)

            for t in range(len(starts3)-1):
                file2 = get_file(f,e,starts3[t+1],ends3[t+1])
                tasmax2_ds = xr.open_dataset(file2)
                tasmax2 = tasmax2_ds['tasmax']
                tasmax2_crop = tasmax2.where(mask_lon & mask_lat, drop=True)
                tasmax1_crop = xr.concat([tasmax1_crop,tasmax2_crop],'time')
            tasmax1_crop.to_netcdf(to_file_path)
            print(f,e)
            






