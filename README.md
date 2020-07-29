# WFCageStudy
Experimental cages were deployed in Shinnecock Bay during the summer of 2016 and 2017.  We monitored growth and survival of YOY winter flounder and simultaneous water quality. Code here is for survival analysis and linear mixed effects models on growth. 

#survivalrevisited2019.R#
- univariate plots
- schoenfield plots
- survival model selection

#growthrevisited2019.R#
-model selection for linear mixed effects models for growth

#growthgraphsOct2019.R#
- generates the visualizations for the final growth models

#survgraphs10_17_19.R#
- generates the visualizations for the final survival models. 

#surv2017_dataminging.R#
- takes indexed_all17_10_17.csv 
- imputes missing data
- converts percent saturation DO to mg/L DO in some cases
- creates time series of environmental data
- creates the vert2_2017.csv
- creates the vert4_2017.csv
- skewness and kurtosis by week
- requires cagetotals_17fixed.csv
- creates the svrt2017forgrowth.csv, which I have renamed to svrt2017forgrowth_regenerated.csv so that it doesn't overwrite the version we currently have. so svrt2017forgrowth.csv is the older version. 
- creates the "weeklysummarydata_2017.csv" which I have renamed to "weeklysummarydata_2017_regenerated.csv" so as to not overwrite the version we have been using, although they should be identical. 
- creates the scaled dataset "weeklysummScaledata_2017.csv", which we don't need, so I silenced it. 
- creates the cumulative summary dataset "cumsummarydata_2017.csv", , which we don't need, so I silenced it. 
- there was some additional exploratory analysis for survival models, but I deleted that in this github version. 

#surv16_dataminging.R#
- takes "Indexed_envirodata16_v2.csv" 
- impute missing data
- converts percent saturation DO to mg/L DO in some cases
- creates time series of environmental data
- creates vert2_2016.csv
- creates vert4_2016.csv
- creates "svrt2016forgrowth.csv" which i renamed to "svrt2016forgrowth_regenerated.csv" so as not to overwrite the existing copy, though they should be identical. 
- creates the "weeklysummarydata_2016_10-8-19.csv" which I have renamed to "weeklysummarydata_2016_10-8-19_regenerated.csv" so as to not overwrite the old version.
- creates "weeklysummarydata_2016_extraweek.csv" which I have silenced because we are not using it. 
- creates the scaled dataset "weeklysummScaledata_2016.csv", which we don't need, so I silenced it. 
- creates the cumulative summary dataset "cumsummarydata_2016.csv", , which we don't need, so I silenced it. 
- - creates the cumulative summary dataset "cumsum2016_deathweekonly.csv", , which we don't need, so I silenced it. 






