# Python codes that support the manuscript "Anthropogenic changes in atmospheric circulation patterns conducive to summer hot extremes in the Northern Hemisphere"
> by Zeqin Huang, Xuezhi Tan\*, Thian Yew Gan, Bingjun Liu\*, Xiaohong Chen

We analyze thermodynamic and dynamic responses of hot extremes to anthropogenic changes in atmospheric circulation patterns. Changes in atmospheric teleconnection patterns under various forcings are interpreted by the self-organizing maps (a python package, [MiniSom](https://github.com/JustGlowing/minisom) is applied). Detection and attribution analyses are performed through the [regularized optimal fingerprinting and the Ribes' attribution method](https://github.com/rafaelcabreu/attribution). 


## Organization of repository
```
|--HWdna_2022
| |--proc_scripts
| |--data_preproc
| |--README.md
| |--Figures
| |--Notebooks
| |--attribution_data
```

## Subdirectory
### data_preproc
>> Preprocessing codes for the raw datasets

```
|--data_preproc
| |--MIROC6_mx2t_preproc.py
| |--IPSL-CM6A-LR_hgt_preproc_500hpa_piControl.py
| |--MRI-ESM2-0_hgt_trend_1979-2014.py
| |--era5_mx2t_preproc.py
| |--IPSL-CM6A-LR_hgt_trend_1979-2014.py
| |--MRI-ESM2-0_mx2t_trend_1979-2014.py
| |--IPSL-CM6A-LR_mx2t_trend_1979-2014.py
| |--MRI-ESM2-0_hgt_preproc_500hpa.py
| |--HadGEM-GC31-LL_hgt_trend_1979-2014.py
| |--ncep_hgt_preproc_500hpa.py
| |--CanESM5_hgt_preproc_500hpa.py
| |--era5_hgt_preproc_500hpa.py
| |--HadGEM-GC31-LL_hgt_preproc_500hpa.py
| |--CanESM5_hgt_trend_1979-2014.py
| |--ncep_preproc_dealwith_2008.py
| |--ncep_tmax_preproc.py
| |--CanESM5_mx2t_trend_1979-2014.py
| |--HadGEM-GC31-LL_hgt_preproc_500hpa_piControl.py
| |--MIROC6_hgt_preproc_500hpa_piControl.py
| |--CanESM5_hgt_preproc_500hpa_piControl.py
| |--MIROC6_hgt_trend_1979-2014.py
| |--ncep_hgt_preproc_200hpa.py
| |--ACCESS-ESM1-5_hgt_preproc_500hpa.py
| |--IPSL-CM6A-LR_hgt_preproc_500hpa.py
| |--MIROC6_hgt_preproc_500hpa.py
| |--MRI-ESM2-0_hgt_preproc_500hpa_piControl.py
| |--era5_hgt_preproc_200hpa.py
| |--HadGEM-GC31-LL_mx2t_trend_1979-2014.py
| |--MIROC6_mx2t_trend_1979-2014.py
```

### proc_scripts
>> Codes for analyses

```
|--proc_scripts
| |--calculate_hot_extreme_occur_WNA.py
| |--determine_winner_for_forcings_relative_to_reanalyses_mean_SOM.py
| |--hgt_tmax_regional_variation.py
| |--concat_hgt_of_forcings_for_SOM.py
| |--calculate_hot_extreme_occur_all_pattern_yearly_piControl.py
| |--calculate_500hpa_GPH_for_patterns_under_forcings.py
| |--calculate_hgt_yearly_for_DNA_piControl.py
| |--calculate_surface_Tmax_historical_average.py
| |--calculate_hot_extreme_occur_trend_sig_concat.py
| |--calculate_tmax_historical_seasonal_cycle_threshold.py
| |--calculate_hot_extreme_occur_all_pattern_yearly.py
| |--som_winner_pattern_historical_nativegrid.py
| |--calculate_surface_Tmax_for_patterns_under_forcings.py
| |--calculate_hgt_yearly_for_DNA.py
| |--calculate_hot_extreme_occur_EAS.py
| |--hot_extreme_per_pattern_occur_trend_concat.py
| |--determine_winner_for_historical_relative_to_reanalyses_mean_SOM.py
| |--calculate_500hpa_GPH_historical_average.py
| |--calculate_hgt_historical_seasonal_cycle.py
| |--calculate_hot_extreme_occur_reanalyses.py
| |--calculate_hot_extreme_occur_EU.py
```
### Notebooks
>> Jupyter notebooks for analyses and visualization

```
|--Notebooks
| |--Optimal_fingerprinting_GPH.ipynb
| |--Fig1_hgt_trend_and_variation.ipynb
| |--Fig2_target_patterns_and_occurrence_under_forcings.ipynb
| |--Fig3_trends_circulation_pattern_and_hot_extreme.ipynb
| |--Fig4_detection_and_attribution_of_hot_extreme.ipynb
| |--Fig5_future_changes_in_patterns_and_hot_extremes.ipynb
| |--FigS1_best_grid_compare.ipynb
| |--FigS2_Trends_of_JJA_GPH_and_Tmax_for_different_forcings.ipynb
| |--FigS3-5_circulation_patterns_for_individual_reanalysis.ipynb
| |--FigS6-11_trends_of_patterns_and_hot_extreme_for_each_subregions.ipynb
| |--FigS12_circulation_anomalies_under_extreme_pattern.ipynb
| |--FigS13_detection_and_attribution_for_reanalyses.ipynb
| |--FigS14_partitioned_trends_in_hot_extreme.ipynb
```
### Figures
>> Rendered figures

```
|--Figures
| |--Fig1_changes_in_GPH_and_associated_hot_extreme.pdf
| |--Fig2_distribution_patt_occur_external_forcing.pdf
| |--Fig3_changes_in_pattern_occurrence_and_associated_hot_extreme.pdf
| |--Fig4_distribution_patt_occur_external_forcing_GPH.pdf
| |--Fig5_future_changes_in_pattern_occurrence_and_associated_hot_extreme.pdf
| |--FigS1_best_grid_selected_for_SOM_analysis.pdf
| |--FigS2_trends_in_circulation_patterns_and_hot_extreme_under_all_forcings.pdf
| |--FigS3_circulation_patterns_categorization_era5.pdf
| |--FigS4_circulation_patterns_categorization_jra55.pdf
| |--FigS5_circulation_patterns_categorization_ncep2.pdf
| |--FigS6_pattern_trend_EU.pdf
| |--FigS7_hot_extreme_trend_EU.pdf
| |--FigS8_pattern_trend_EAS.pdf
| |--FigS9_hot_extreme_trend_EAS.pdf
| |--FigS10_pattern_trend_WNA.pdf
| |--FigS11_hot_extreme_trend_WNA.pdf
| |--FigS12_circulation_anomalies_under_external_forcings.pdf
| |--FigS13_detection_and_attribution_for_reanalysis_using_ROF_Ribes.pdf
| |--FigS14_partitioned_trends_in_hot_extreme.pdf
```

## Datasets used
- Detection and Attribution Model Intercomparison (DAMIP) dataset from [CMIP6](https://esgf-node.llnl.gov/search/cmip6/)
- [ERA5](https://cds.climate.copernicus.eu) from European Centre for Medium-Range Weather Forecasts (ECMWF)
- [NCEP2](https://www.esrl.noaa.gov/psd/data/gridded/) provided by the NOAA/OAR/ESRL PSL
- [JRA55](https://climatedataguide.ucar.edu/climate-data/jra-55) available at the National Center for Atmospheric Research

