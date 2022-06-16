'''
This script is for concatenating the extracted 'hot extreme per pattern occurrence' for different models under different external forcings
'''
import pandas as pd

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

def get_per_pattern_occur(dataset_name, forcing, run, domain):
    path = '/home/xtan/scratch/hzq/HWdna/procData/' + dataset_name + '_' + forcing + '_' + run + '_hot_extreme_per_pattern_occur_' + domain + '.csv'
    per_pattern_occur = pd.read_csv(path,index_col=0)
    return per_pattern_occur

def get_per_pattern_occur_reanalyses(dataset_name, domain):
    path = '/home/xtan/scratch/hzq/HWdna/procData/' + dataset_name + '_hot_extreme_per_pattern_occur_' + domain + '.csv'
    per_pattern_occur = pd.read_csv(path,index_col=0)
    return per_pattern_occur

def save_results(forcing,domain,patt):## patt must be one of [0,1,2,3]
    patt_name = ['Pattern1','Pattern2','Pattern3','Pattern4']
    per_patter_occur_cmip_df = get_per_pattern_occur(dataset_name='CanESM5',forcing=forcing,run='r1i1p1f1',domain=domain)
    per_patter_occur_cmip_df = per_patter_occur_cmip_df[patt_name[patt]]
    per_patter_occur_cmip_df.name = 'CanESM5_r1i1p1f1'
    for d in dataset_src_run.keys():
        for r in dataset_src_run[d][forcing]:
            if d == 'CanESM5' and r == 'r1i1p1f1':
                continue
            else:
                per_pattern_occur = get_per_pattern_occur(dataset_name=d,forcing=forcing,run=r,domain=domain)
                per_pattern_occur = per_pattern_occur[patt_name[patt]]
                per_pattern_occur.name = d + '_' + r
                per_patter_occur_cmip_df = pd.concat([per_patter_occur_cmip_df,per_pattern_occur],axis=1)
                print(d + '_' + r)
    to_path = '/home/xtan/scratch/hzq/HWdna/procData/' + 'hot_extreme_per_pattern_occur_patt' + str(patt + 1) + '_variation_' + forcing + '_' + domain + '.csv'
    per_patter_occur_cmip_df.to_csv(to_path)

def save_results_reanalyses(domain,patt):## patt must be one of [0,1,2,3]
    patt_name = ['Pattern1','Pattern2','Pattern3','Pattern4']
    per_patter_occur_cmip_df = get_per_pattern_occur_reanalyses(dataset_name='era5',domain=domain)
    per_patter_occur_cmip_df = per_patter_occur_cmip_df[patt_name[patt]]
    per_patter_occur_cmip_df.name = 'era5'
    for d in ['era5','jra55','ncep2']:
        if d == 'era5':
            continue
        else:
            per_pattern_occur = get_per_pattern_occur_reanalyses(dataset_name=d,domain=domain)
            per_pattern_occur = per_pattern_occur[patt_name[patt]]
            per_pattern_occur.name = d
            per_patter_occur_cmip_df = pd.concat([per_patter_occur_cmip_df,per_pattern_occur],axis=1)
            print(d)
    to_path = '/home/xtan/scratch/hzq/HWdna/procData/' + 'hot_extreme_per_pattern_occur_patt' + str(patt + 1) + '_variation_' + 'reanalyses' + '_' + domain + '.csv'
    per_patter_occur_cmip_df.to_csv(to_path)

for f in ['historical','hist-GHG','hist-aer','hist-nat','ssp585']:
    for d in ['EU','EAS','WNA']:
        for p in [0,1,2,3]:
            save_results(forcing=f,domain=d,patt=p)

for d in ['EU','EAS','WNA']:
    for p in [0,1,2,3]:
        save_results_reanalyses(domain=d,patt=p)
