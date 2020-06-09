# Franzen_InfoFlows_2020
python codes for Franzen, Farahani, Goodwell et al, WRR 2020 (in review)

DOI for this repository: 10.5281/zenodo.3860555
https://zenodo.org/badge/latestdoi/267372638

Main analysis codes:
COHW_analysis_overall.py computes IT measures for multi-year time window between each precipitation grid cell and outlet flow rates, saves as excel file

COHW_analysis_annual_trends.py computes IT measures for annual windows for each season, detects trends using linear regression and Sen slope methods, saves trend results in excel files


Functions: compute_info_measures.py and compute_icrit.py for information theory calls

Datasets: 
CO_rainfall_data.pickle is a subset of CPC Gage-based gridded precipitation data for Colorado Headwaters HUC4 basin (CPC US Unified Precipitation data provided by the NOAA/OAR/ESRL PSL, Boulder, Colorado, USA, from their Web site at https://psl.noaa.gov/)

COHW Gage Data 1952-2017.csv is from USGS Water Data (https://waterdata.usgs.gov/nwis) for daily flow data, Colorado River near Utah State Line, and Gunnison River near Grand Junction gages
