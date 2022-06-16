import xarray as xr

# ds1 = xr.open_dataset('vwnd.10m.gauss.2008.nc')
ds2 = xr.open_dataset('uwnd.10m.gauss.2008.nc')
ds3 = xr.open_dataset('air.2m.gauss.2008.nc')
ds4 = xr.open_dataset('shum.2m.gauss.2008.nc')
ds5 = xr.open_dataset('tcdc.eatm.gauss.2008.nc')
ds6 = xr.open_dataset('prate.sfc.gauss.2008.nc')

# ds1 = ds1.resample(time='D').mean()
ds2 = ds2.resample(time='D').mean()
ds3 = ds3.resample(time='D').mean()
ds4 = ds4.resample(time='D').mean()
ds5 = ds5.resample(time='D').mean()
ds6 = ds6.resample(time='D').mean()

# ds1.to_netcdf('vwnd.10m.gauss.2008.nc')
ds2.to_netcdf('uwnd.10m.gauss.2008.nc')
ds3.to_netcdf('air.2m.gauss.2008.nc')
ds4.to_netcdf('shum.2m.gauss.2008.nc')
ds5.to_netcdf('tcdc.eatm.gauss.2008.nc')
ds6.to_netcdf('prate.sfc.gauss.2008.nc')
