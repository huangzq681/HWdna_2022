import xarray as xr
import pandas as pd
import numpy as np
from scipy.stats import linregress

time_start = '1979-01-01'
time_end   = '2014-12-31'

## select level and date range
dataset_src_run = {
    'CanESM5':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1','r8i1p1f1','r9i1p1f1','r10i1p1f1'],
    },
    'HadGEM3-GC31-LL':{
        'historical':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-GHG':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-nat':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'hist-aer':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
        'ssp585':['r1i1p1f3','r2i1p1f3','r3i1p1f3','r4i1p1f3'],
    },
    'MIROC6':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1'],
    },
    'IPSL-CM6A-LR':{
        'historical':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r5i1p1f1','r6i1p1f1','r7i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r4i1p1f1','r6i1p1f1'], #'r3i1p1f1','r14i1p1f1' are excluded, as there are not corresponding runs in historical forcing
    },
    'MRI-ESM2-0':{
        'historical':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-GHG':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-nat':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'hist-aer':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
        'ssp585':['r1i1p1f1','r2i1p1f1','r3i1p1f1','r4i1p1f1','r5i1p1f1'],
    },
}

def _compute_slope(var):
    slp = linregress(range(len(var)),var).slope
    return slp

def _compute_sig(var):
    sig = linregress(range(len(var)),var).pvalue
    return sig

def trend_cal(data):
    slopes = xr.apply_ufunc(_compute_slope,
                            data,
                            vectorize=True,
                            dask='parallelized', 
                            input_core_dims=[['year']],
                            output_dtypes=[float],
                            )
    sig = xr.apply_ufunc(_compute_sig,
                            data,
                            vectorize=True,
                            dask='parallelized', 
                            input_core_dims=[['year']],
                            output_dtypes=[float],
                            )
    return slopes,sig

def get_file(dataset_name,domain,pattern,forcing=None,run=None): ## pattern must be one of [0,1,2,3]
    if dataset_name in ['era5','jra55','ncep2']:
        hot_extr_patt_occur_yearly_path = '/lustre07/scratch/xtan/hzq/HWdna/procData/' + dataset_name + '_hot_extreme_patt' + str(pattern + 1) + '_occur_' + domain + '.nc'
    else:
        hot_extr_patt_occur_yearly_path = '/lustre07/scratch/xtan/hzq/HWdna/procData/' + dataset_name + '_' + forcing + '_' + run + '_hot_extreme_patt' + str(pattern + 1) + '_occur_' + domain + '.nc'
    return hot_extr_patt_occur_yearly_path

def to_ncfile(domain,pattern,forcing,trendOrSig = 'trend'):
    to_file_tmax_ano_forcing_trendOrSig = '/lustre07/scratch/xtan/hzq/HWdna/procData/' + 'hot_extreme_patt'+ str(pattern + 1) +'_occur_' + trendOrSig + '_' + forcing + '_' + domain + '.nc'
    return to_file_tmax_ano_forcing_trendOrSig

def to_csvfile(domain,pattern,forcing):
    to_file_tmax_ano_forcing_variation = '/lustre07/scratch/xtan/hzq/HWdna/procData/' + 'hot_extreme_patt'+ str(pattern + 1) +'_occur_variation_' + forcing + '_' + domain + '.csv'
    return to_file_tmax_ano_forcing_variation

def save_results(domain,patt):
    for f in dataset_src_run['CanESM5'].keys():
        if f == 'historical':
            ## spatial trend and significant
            era5_path = get_file(dataset_name='era5',domain=domain,pattern=patt)
            jra55_path = get_file(dataset_name='jra55',domain=domain,pattern=patt)
            ncep2_path = get_file(dataset_name='ncep2',domain=domain,pattern=patt)
            tmax_ano_era5 = xr.open_dataset(era5_path)
            tmax_ano_era5 = tmax_ano_era5.rename({'mx2t':'tasmax'})
            tmax_ano_era5 = tmax_ano_era5['tasmax']
            tmax_ano_jra55 = xr.open_dataset(jra55_path)['tasmax']
            tmax_ano_ncep2 = xr.open_dataset(ncep2_path)
            tmax_ano_ncep2 = tmax_ano_ncep2.rename({'tmax':'tasmax'})
            tmax_ano_ncep2 = tmax_ano_ncep2['tasmax']
            tmax_ano_era5_trend, tmax_ano_era5_sig = trend_cal(tmax_ano_era5)
            tmax_ano_jra55_trend, tmax_ano_jra55_sig = trend_cal(tmax_ano_jra55)
            tmax_ano_ncep2_trend, tmax_ano_ncep2_sig = trend_cal(tmax_ano_ncep2)
            tmax_ano_era5_trend = tmax_ano_era5_trend.expand_dims({'Dataset_run':1})
            tmax_ano_jra55_trend = tmax_ano_jra55_trend.expand_dims({'Dataset_run':1})
            tmax_ano_ncep2_trend = tmax_ano_ncep2_trend.expand_dims({'Dataset_run':1})
            tmax_ano_era5_sig = tmax_ano_era5_sig.expand_dims({'Dataset_run':1})
            tmax_ano_jra55_sig = tmax_ano_jra55_sig.expand_dims({'Dataset_run':1})
            tmax_ano_ncep2_sig = tmax_ano_ncep2_sig.expand_dims({'Dataset_run':1})
            tmax_ano_era5_trend = tmax_ano_era5_trend.assign_coords(Dataset_run = ['era5'])
            tmax_ano_jra55_trend = tmax_ano_jra55_trend.assign_coords(Dataset_run = ['jra55'])
            tmax_ano_ncep2_trend = tmax_ano_ncep2_trend.assign_coords(Dataset_run = ['ncep2'])
            tmax_ano_era5_sig = tmax_ano_era5_sig.assign_coords(Dataset_run = ['era5'])
            tmax_ano_jra55_sig = tmax_ano_jra55_sig.assign_coords(Dataset_run = ['jra55'])
            tmax_ano_ncep2_sig = tmax_ano_ncep2_sig.assign_coords(Dataset_run = ['ncep2'])
            tmax_ano_trend_all = xr.concat([tmax_ano_era5_trend,tmax_ano_jra55_trend,tmax_ano_ncep2_trend],dim='Dataset_run')
            tmax_ano_sig_all = xr.concat([tmax_ano_era5_sig,tmax_ano_jra55_sig,tmax_ano_ncep2_sig],dim='Dataset_run')
            ## spatial weighted mean temporal variation
            weights = np.cos(np.deg2rad(tmax_ano_era5.lat))
            weights.name = "weights"
            tmax_ano_era5_weight = tmax_ano_era5.weighted(weights)
            tmax_ano_era5_weight = tmax_ano_era5_weight.mean(("lon", "lat"))
            tmax_ano_era5_weight = tmax_ano_era5_weight.to_series()
            tmax_ano_era5_weight.name = 'era5'
            tmax_ano_jra55_weight = tmax_ano_jra55.weighted(weights)
            tmax_ano_jra55_weight = tmax_ano_jra55_weight.mean(("lon", "lat"))
            tmax_ano_jra55_weight = tmax_ano_jra55_weight.to_series()
            tmax_ano_jra55_weight.name = 'jra55'
            tmax_ano_ncep2_weight = tmax_ano_ncep2.weighted(weights)
            tmax_ano_ncep2_weight = tmax_ano_ncep2_weight.mean(("lon", "lat"))
            tmax_ano_ncep2_weight = tmax_ano_ncep2_weight.to_series()
            tmax_ano_ncep2_weight.name = 'ncep2'
            # tmax_ano_variation_all = pd.DataFrame({'era5':tmax_ano_era5_weight.values,'jra55':tmax_ano_jra55_weight.values,'ncep2':tmax_ano_ncep2_weight})
            tmax_ano_variation_all = pd.concat([tmax_ano_era5_weight,tmax_ano_jra55_weight,tmax_ano_ncep2_weight], axis=1)

            for d in dataset_src_run.keys():
                for r in dataset_src_run[d][f]:
                    ## spatial trend and significant
                    dr_path = get_file(dataset_name=d,pattern=patt,domain=domain,forcing=f,run=r)
                    tmax_ano = xr.open_dataset(dr_path)['tasmax']
                    tmax_ano_trend, tmax_ano_sig = trend_cal(tmax_ano)
                    tmax_ano_trend = tmax_ano_trend.expand_dims({'Dataset_run':1})
                    tmax_ano_trend = tmax_ano_trend.assign_coords(Dataset_run = [d + '_' + r])
                    tmax_ano_sig = tmax_ano_sig.expand_dims({'Dataset_run':1})
                    tmax_ano_sig = tmax_ano_sig.assign_coords(Dataset_run = [d + '_' + r])
                    tmax_ano_trend_all= xr.concat([tmax_ano_trend_all,tmax_ano_trend],dim='Dataset_run')
                    tmax_ano_sig_all= xr.concat([tmax_ano_sig_all,tmax_ano_sig],dim='Dataset_run')
                    ## spatial weighted mean temporal variation
                    tmax_ano_weight = tmax_ano.weighted(weights)
                    tmax_ano_weight = tmax_ano_weight.mean(("lon", "lat"))
                    tmax_ano_weight = tmax_ano_weight.to_series()
                    tmax_ano_weight.name = d + '_' + r
                    tmax_ano_variation_all = pd.concat([tmax_ano_variation_all,tmax_ano_weight], axis=1)
                    print(f+'_'+d+'_'+r)
            tmax_ano_trend_all.to_netcdf(to_ncfile(domain=domain,pattern=patt,forcing=f,trendOrSig='trend'))
            tmax_ano_sig_all.to_netcdf(to_ncfile(domain=domain,pattern=patt,forcing=f,trendOrSig='sig'))
            tmax_ano_variation_all.to_csv(to_csvfile(domain=domain,pattern=patt,forcing=f))

        else:
            ## spatial trend and significant
            d = 'CanESM5'
            r = dataset_src_run[d][f][0]
            dr_path = get_file(dataset_name=d,pattern=patt,domain=domain,forcing=f,run=r)
            tmax_ano = xr.open_dataset(dr_path)['tasmax']
            tmax_ano_trend, tmax_ano_sig = trend_cal(tmax_ano)
            tmax_ano_trend = tmax_ano_trend.expand_dims({'Dataset_run':1})
            tmax_ano_trend = tmax_ano_trend.assign_coords(Dataset_run = [d + '_' + r])
            tmax_ano_sig = tmax_ano_sig.expand_dims({'Dataset_run':1})
            tmax_ano_sig = tmax_ano_sig.assign_coords(Dataset_run = [d + '_' + r])
            tmax_ano_trend_all= tmax_ano_trend
            tmax_ano_sig_all= tmax_ano_sig
            ## spatial weighted mean temporal variation
            weights = np.cos(np.deg2rad(tmax_ano.lat))
            weights.name = "weights"
            tmax_ano_weight = tmax_ano.weighted(weights)
            tmax_ano_weight = tmax_ano_weight.mean(("lon", "lat"))
            tmax_ano_variation_all = tmax_ano_weight.to_series()

            for d in dataset_src_run.keys():
                for r in dataset_src_run[d][f]:
                    ## spatial trend and significant
                    if d == 'CanESM5' and r == dataset_src_run[d][f][0]:
                        continue
                    else:
                        dr_path = get_file(dataset_name=d,pattern=patt,domain=domain,forcing=f,run=r)
                        tmax_ano = xr.open_dataset(dr_path)['tasmax']
                        tmax_ano_trend, tmax_ano_sig = trend_cal(tmax_ano)
                        tmax_ano_trend = tmax_ano_trend.expand_dims({'Dataset_run':1})
                        tmax_ano_trend = tmax_ano_trend.assign_coords(Dataset_run = [d + '_' + r])
                        tmax_ano_sig = tmax_ano_sig.expand_dims({'Dataset_run':1})
                        tmax_ano_sig = tmax_ano_sig.assign_coords(Dataset_run = [d + '_' + r])
                        tmax_ano_trend_all= xr.concat([tmax_ano_trend_all,tmax_ano_trend],dim='Dataset_run')
                        tmax_ano_sig_all= xr.concat([tmax_ano_sig_all,tmax_ano_sig],dim='Dataset_run')
                        ## spatial weighted mean temporal variation
                        tmax_ano_weight = tmax_ano.weighted(weights)
                        tmax_ano_weight = tmax_ano_weight.mean(("lon", "lat"))
                        tmax_ano_weight = tmax_ano_weight.to_series()
                        tmax_ano_weight.name = d + '_' + r
                        tmax_ano_variation_all = pd.concat([tmax_ano_variation_all,tmax_ano_weight], axis=1)
                        print(f+'_'+d+'_'+r) 
            tmax_ano_trend_all.to_netcdf(to_ncfile(domain=domain,pattern=patt,forcing=f,trendOrSig='trend'))
            tmax_ano_sig_all.to_netcdf(to_ncfile(domain=domain,pattern=patt,forcing=f,trendOrSig='sig'))
            tmax_ano_variation_all.to_csv(to_csvfile(domain=domain,pattern=patt,forcing=f))       

save_results(domain='EU',patt=0)
save_results(domain='EU',patt=1)
save_results(domain='EU',patt=2)
save_results(domain='EU',patt=3)

save_results(domain='EAS',patt=0)
save_results(domain='EAS',patt=1)
save_results(domain='EAS',patt=2)
save_results(domain='EAS',patt=3)

save_results(domain='WNA',patt=0)
save_results(domain='WNA',patt=1)
save_results(domain='WNA',patt=2)
save_results(domain='WNA',patt=3)