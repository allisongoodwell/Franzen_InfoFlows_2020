# -*- coding: utf-8 -*-
"""
Updated 6/5/2020

This is modified from Sam Franzen's original code, to input flow and precip data
first compute information measures year by year, select 1 year windows, or 3 year windows
then compute trends in different meausures: mutual information, entropy, mutual information divided by streamflow entropy

Trends are computed for measures of normalized mutual information (I/H(Q)), lag time (L), threshold value (T)

Inputs: CPC gage-based gridded data for CO Headwaters HUC4 basin region, USGS flow rate data

#functions: comput_info_measures.py and compute_icrit.py for IT computations

Outputs: excel files with trend values (statistically signficant trends) for both linear regression and Sen Slope method

@author: Sam, Allison, Mozhgan
"""


import pickle
import numpy as np
import pandas as pd
import time
from itertools import chain
from compute_icrit import compute_icrit as ci
from compute_info_measures import compute_info_measures as cim
from scipy import stats
import matplotlib.pyplot as plt

  

start = time.time()
 
start_year1 = 1953
end_year1 = 2018


#designate bins for P (precip), xbin, and Q (flow rate), ybin
xbin = 2
ybin = 5

#designate window of years (should be odd) for each year analysis
window = 1

#designate max lag time (days)
delay = 7

#minimum threshold values (mm) for daily precipitation
thresholds = [0.3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]


#create year, month, day vectors
ndays_per_month = np.asfarray([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

yearlist = []
monthlist = []
daylist = []
seasonlist = []

#define seasons
for y in range (start_year1-1,end_year1):
    for m in range (0,12):
        if y in range (1952, 2018, 4) and m == 1 or y == 1952 and m == 1:
            length  = int(29)
        else:
            length = int(ndays_per_month[m])
        
        if m in range(0,2) or m ==11: #winter
            s = 1
        elif m in range(2,5):
            s = 2
        elif m in range(5,8):
            s = 3
        else:
            s = 4
            
        for d in range (1,length+1):
            daylist.append(d)
            monthlist.append(m+1)
            yearlist.append(y)
            seasonlist.append(s)
            
dfyear = pd.DataFrame(yearlist, columns = ['year'])  
dfmonth = pd.DataFrame(monthlist, columns = ['month']) 
dfday = pd.DataFrame(daylist, columns = ['day']) 
dfseason = pd.DataFrame(seasonlist, columns = ['season']) 

dfdate = pd.concat([dfyear, dfmonth, dfday, dfseason], axis = 1)

dateyrarray = np.array([[yearlist]])

#load daily streamflow data (difference between COHW and Gunnison flow rates)
dfsf1 = pd.read_csv('COHW Gage Data 1952-2017.csv', usecols=[4])
dfsf1.rename(columns={'difference': 'flow'}, inplace=True)
dfsf1_list = dfsf1['flow']
dfsf1_list = dfsf1_list.tolist()
dfsf1_series = pd.Series(dfsf1_list)
dfsf1_array = np.array(dfsf1_series)

#load gridded daily rainfall data
f_myfile = open('CO_rainfall_data.pickle','rb')
Co_pptdata = pickle.load(f_myfile) #ppt_long_list = rainfall data [0] and percentiles [1]
lat = pickle.load(f_myfile)
lon = pickle.load(f_myfile)
f_myfile.close()

dfppt = pd.DataFrame(Co_pptdata)




#%%


yearvect = range(start_year1,end_year1,1)

#create lists and arrays for data storage
winterplist = []
springplist = []
summerplist = []
fallplist = []

winterflist = []
springflist = []
summerflist = []
fallflist = []

winterI = []
springI = []
summerI = []
fallI = []

wintersgnf = []
springsgnf = []
summersgnf = []
fallsgnf = []

winterlen = np.sum([31,31,28]) #DJF
springlen = np.sum([31,30,31]) #MAM
summerlen = np.sum([30, 31, 31]) #JJA
falllen = np.sum([30, 31, 30]) #SON

#breakpoints between years within a 5-year period
winter_breakpoints = [winterlen*i for i in range(1,window)]
spring_breakpoints = [springlen*i for i in range(1,window)]
summer_breakpoints = [summerlen*i for i in range(1,window)]
fall_breakpoints = [falllen*i for i in range(1,window)]


winterIp = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
springIp = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
summerIp = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
fallIp = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))

winterHQ = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
springHQ =  np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
summerHQ =  np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
fallHQ =  np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))

winterSI = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
springSI = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
summerSI = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
fallSI = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))

winterp1 = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
springp1 = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
summerp1 = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))
fallp1 = np.zeros((len(Co_pptdata),len(yearvect),len(thresholds),delay+1))


#load ppt data, specific station
#for each gridded cell...
for pg in range (0, len(Co_pptdata)):
    pptr8 = []
    
    print(pg)
    
    pptr8 = list(chain(*Co_pptdata[pg]))
    del pptr8[0:1461]
    dfpptr8 = pd.DataFrame(pptr8)

    dfppsf = pd.concat([dfpptr8, dfsf1], axis = 1)
      
    #for each threshold magnitude value....
    for thvect,th in enumerate(thresholds):
        
        plimit = th
        
        bi_ppt_list = []
    
        #check for NaN values in each data point....
        for i in range(len(dfppsf)):
                        
            p1 = dfppsf.iloc[i, 0]
            if p1 == 'NaN':
                p1 = 'NaN'
            elif p1 > plimit:
                p1 = 1
            elif p1 <= plimit:
                p1 = 0
            bi_ppt_list.append(p1)
        
        dfpptr8 = pd.DataFrame(bi_ppt_list, columns=['ppt8'])    
        
        dfcc = pd.concat([dfdate, dfpptr8, dfsf1], axis = 1)
        dfcc = dfcc.dropna(0)
        dfcc = dfcc.set_index(['year'])                                
        dfcc['year'] = dfcc.index
        
        di_list = []
        
        pptlist = bi_ppt_list
        sflist = dfppsf['flow'].tolist()
    
#AEG: I'm changing this loop to go by ones instead of 3 -- to get 3 year "moving window" instead of one value every 5 years
#now y_ind is the index (year 0, 1, 2, ...), g is the actual year (1950, 1951....))  
        #for each year....
        for y_ind,g in enumerate(yearvect):
            start_yearmid = g
            start_year = int(g-((window-1)/2))
            
            end_year = start_year + window
            
            #AEG: added this so that the plist and flists don't keep expanding on top of eachother (as it was, they get longer and longer)
            springplist=[]
            springflist=[]
            fallplist=[]
            fallflist=[]
            summerplist=[]
            summerflist=[]
            winterplist=[]
            winterflist=[]
        
            
            for y in range(start_year, end_year):
                
                for s in range(1,5):
                    if s>1:  #spring, summer, fall get data in X year window
                        df_short = dfcc[np.logical_and(dfcc['season']==s,dfcc['year']==y)]
                        plist = df_short['ppt8'].tolist()
                        flist = df_short['flow'].tolist()
        
                        #AEG: something I realized - after this loop, we are appending (extending) the list to get string of data that is X years
                        #worth of data points.  Due to season splitting, some of these data points are not actually continuous.
                        #So when we do lagged pdfs, some of the points are "cheating" - lagging from one season to the next...
                        #For 0 lag, it is fine, but for any other lag, need to leave out certain points...
                        
                        if s == 2:
                            springplist.extend(plist)
                            springflist.extend(flist)
            
                        if s == 3:
                            summerplist.extend(plist)
                            summerflist.extend(flist)
                            
                        if s == 4:
                            fallplist.extend(plist)
                            fallflist.extend(flist)
                        
                    else: #winter is awkward (want Dec of past year and Jan-Feb of this year) get data in 3 year window
                            #data_indices = [i for i in range(0,len(dfcc)) if np.logical_or((np.logical_and(dfcc.loc[i]['month'] in range(1,3), dfcc.loc[i]['year']==y)), (np.logical_and(dfcc.loc[i]['month']==12, dfcc.loc[i]['year']==y-1)))]
                            
                        df_short = dfcc[np.logical_or(np.logical_and(dfcc['month']<3,dfcc['year']==y),np.logical_and(dfcc['month']==12,dfcc['year']==y-1))]
                        plist = df_short['ppt8'].tolist()
                        flist = df_short['flow'].tolist()
                        
                        winterplist.extend(plist)
                        winterflist.extend(flist)
            
            #loop through seasons
            for t in range (1,5):
            
                if t == 1:
                    pptlist = winterplist
                    sflist = winterflist
                    breakpoints = winter_breakpoints
                    
                if t == 2:
                    pptlist = springplist
                    sflist = springflist
                    breakpoints = spring_breakpoints
                    
                if t == 3:
                    pptlist = summerplist
                    sflist = summerflist
                    breakpoints = summer_breakpoints
                    
                if t == 4:
                    pptlist = fallplist
                    sflist = fallflist
                    breakpoints = fall_breakpoints
                    
                #loop through lag times
                for d in range(0, delay+1):
                    
                    pptlist_new = pptlist[0:len(pptlist)-d] 
                    sflist_new = sflist[d:len(sflist)]
                    
                    #AEG: here is the strange part where we have to leave out certain points...(breakpoints)
                    if d>0:
                        for b in sorted(breakpoints,reverse=True):
                            pptlist_new[(b-d+1):b]=[]
                            sflist_new[(b-d+1):b]=[]
                    
                    hist_array = np.array((pptlist_new, sflist_new)).T
                    pdf, edges = np.histogramdd(hist_array, (xbin, ybin))
                    pdf = pdf/np.sum(pdf)
                    I_dict = cim(pdf)
                    Icrit_dict = ci(pptlist_new, sflist_new)
                    
                    if t == 1:
                        a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                        if a > 0:

                            winterIp[pg,y_ind,thvect,d]=I_dict['Ix1x2']
                            winterSI[pg,y_ind,thvect,d]=I_dict['SI']
                            winterp1[pg,y_ind,thvect,d]=I_dict['m1']                        
                            winterHQ[pg,y_ind,thvect,d]=I_dict['Hx2']
                                        
                    if t == 2:
                        a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                        if a > 0:

                            springIp[pg,y_ind,thvect,d]=I_dict['Ix1x2']
                            springSI[pg,y_ind,thvect,d]=I_dict['SI']
                            springp1[pg,y_ind,thvect,d]=I_dict['m1']                       
                            springHQ[pg,y_ind,thvect,d]=I_dict['Hx2']

                    
                    if t == 3:
                        a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                        if a > 0:

                            summerIp[pg,y_ind,thvect,d]=I_dict['Ix1x2']
                            summerSI[pg,y_ind,thvect,d]=I_dict['SI']
                            summerp1[pg,y_ind,thvect,d]=I_dict['m1']                        
                            summerHQ[pg,y_ind,thvect,d]=I_dict['Hx2']
                    
                    if t == 4:
                        a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                        if a > 0:

                            fallIp[pg,y_ind,thvect,d]=I_dict['Ix1x2']
                            fallSI[pg,y_ind,thvect,d]=I_dict['SI']
                            fallp1[pg,y_ind,thvect,d]=I_dict['m1']                       
                            fallHQ[pg,y_ind,thvect,d]=I_dict['Hx2']


end = time.time()
print ((end-start)/60)
print('done with main IT analysis, now detecting trends...')

#%%
#next, find trends in these annual values
#create arrays and lists for data storage


maxI = np.zeros((len(Co_pptdata),4,len(yearvect)))
maxL = np.zeros((len(Co_pptdata),4,len(yearvect)))
maxT = np.zeros((len(Co_pptdata),4,len(yearvect)))
maxHQ = np.zeros((len(Co_pptdata),4,len(yearvect)))

SenSlope_I=np.zeros((len(Co_pptdata),4))
SenSlope_L=np.zeros((len(Co_pptdata),4))
SenSlope_T=np.zeros((len(Co_pptdata),4))
SenSlope_H=np.zeros((len(Co_pptdata),4))

LinRegSlope_I=np.zeros((len(Co_pptdata),4))
LinRegSlope_L=np.zeros((len(Co_pptdata),4))
LinRegSlope_T=np.zeros((len(Co_pptdata),4))
LinRegSlope_H=np.zeros((len(Co_pptdata),4))



singleyearlist = list(yearvect)
x = np.asarray(singleyearlist)

CI=.95 #statistical significance level

#conduct data analysis to find max mutual info and corresponding threshold and lag for each year in each season
#find trend(slope) of I values across 60 year period
for pg in range (0, len(Co_pptdata)):

    for season in range(0,4):
        
        if (season==0):
            Ip = winterIp
            HQ = winterHQ
        elif (season==1):
            Ip = springIp
            HQ = springHQ
        elif (season==2):
            Ip = summerIp
            HQ = summerHQ
        elif (season==3):
            Ip = fallIp
            HQ = fallHQ
            
        
        #iterate through year windows, find dominant information values, lags, thresholds
        for y_ind,g in enumerate(yearvect):
            
            maxI[pg,season,y_ind]= np.max(Ip[pg,y_ind].flatten(), axis=0)
            maxthresh = int(np.ceil(np.argmax(Ip[pg,y_ind])/8))
            if maxthresh == 0:
                maxthresh = 0
            else: 
                maxthresh = maxthresh-1
            if maxthresh == 0:
                maxthresh1 = 0.3
            else: 
                maxthresh1 = maxthresh
                
            maxT[pg,season, y_ind]=maxthresh1
            maxlag = int(np.argmax(Ip[pg,y_ind,maxthresh]))
            maxL[pg,season,y_ind] = maxlag
            
            maxHQ[pg,season,y_ind] = HQ[pg, y_ind,maxthresh,maxlag]
          
 
        #trend in mutual information
        y = maxI[pg,season]
        
        res = stats.theilslopes(y, x, CI)
        lsq_res, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
        if (res[2] > 0 and res [3] > 0) or (res[2] < 0 and res [3] < 0):
            SenSlope_I[pg,season]=res[0]

        if (p_value < 1-CI):
            LinRegSlope_I[pg,season]=lsq_res

        #trend in lag
        y = maxL[pg,season]
        
        res = stats.theilslopes(y, x, CI)
        lsq_res, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
        if (res[2] > 0 and res [3] > 0) or (res[2] < 0 and res [3] < 0):
            SenSlope_L[pg,season]=res[0]

        if (p_value < 1-CI):
            LinRegSlope_L[pg,season]=lsq_res
            
        #trend in threshold
        y = maxT[pg,season]
        
        res = stats.theilslopes(y, x, CI)
        lsq_res, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
        if (res[2] > 0 and res [3] > 0) or (res[2] < 0 and res [3] < 0):
            SenSlope_T[pg,season]=res[0]

        if (p_value < 1-CI):
            LinRegSlope_T[pg,season]=lsq_res
            
            
        #trend in I/H
        y = np.divide(maxI[pg,season],maxHQ[pg,season])
        
        res = stats.theilslopes(y, x, CI)
        lsq_res, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    
        if (res[2] > 0 and res [3] > 0) or (res[2] < 0 and res [3] < 0):
            SenSlope_H[pg,season]=res[0]

        if (p_value < 1-CI):
            LinRegSlope_H[pg,season]=lsq_res 
           
            
            
        #make some plots for cases with trends   
        if LinRegSlope_I[pg,season]!=0:
            plt.figure(1)
            plt.title('Mutual Information Trends')
            plt.subplot(4,1,season+1)
            plt.plot(maxI[pg,season])
            
        if LinRegSlope_T[pg,season]!=0:
            plt.figure(2)
            plt.title('Threshold Trends')
            plt.subplot(4,1,season+1)
            plt.plot(maxT[pg,season])

        if LinRegSlope_L[pg,season]!=0:
            plt.figure(3)
            plt.title('Lag Trends')
            plt.subplot(4,1,season+1)
            plt.plot(maxL[pg,season])

     
        if LinRegSlope_H[pg,season]!=0:
            plt.figure(4)
            plt.title('I/H Trends')
            plt.subplot(4,1,season+1)
            plt.plot(np.divide(maxI[pg,season],maxHQ[pg,season])*100)  
            
            
#%% save results into a dataframe and excel file

DF_IHtrends = pd.DataFrame(LinRegSlope_H)
DF_Lagtrends = pd.DataFrame(LinRegSlope_L)
DF_Threshtrends = pd.DataFrame(LinRegSlope_T)

DF_IHtrends.to_excel("LinRegTrends_IH_oneyrwindow.xlsx")
DF_Lagtrends.to_excel("LinRegTrends_Lag_oneyrwindow.xlsx")
DF_Threshtrends.to_excel("LinRegTrends_Thresh_oneyrwindow.xlsx")

DF_IHtrends_S = pd.DataFrame(SenSlope_H)
DF_Lagtrends_S = pd.DataFrame(SenSlope_L)
DF_Threshtrends_S = pd.DataFrame(SenSlope_T)

DF_IHtrends_S.to_excel("SenSlope_IH_oneyrwindow.xlsx")
DF_Lagtrends_S.to_excel("SenSlope_Lag_oneyrwindow.xlsx")
DF_Threshtrends_S.to_excel("SenSlope_Thresh_oneyrwindow.xlsx")

#%% find annual average values for each measure (I/H, L, and T), and percent changes

from numpy import inf

maxIoverH = np.true_divide(maxI,maxHQ)
maxIoverH[maxIoverH==inf]=0

avg_I_vals = np.average(maxIoverH,axis=2)
avg_L_vals = np.average(maxL,axis=2)
avg_T_vals = np.average(maxT,axis=2)

avg_I_change = ((LinRegSlope_H*65)/avg_I_vals)*100
avg_L_change = ((LinRegSlope_L*65)/avg_L_vals)*100
avg_T_change = ((LinRegSlope_T*65)/avg_T_vals)*100

avg_percent_IH=np.true_divide(avg_I_change.sum(0),(avg_I_change!=0).sum(0))
avg_percent_Lag=np.true_divide(avg_L_change.sum(0),(avg_L_change!=0).sum(0))

avg_I_vals[avg_I_change==0]=np.nan
avg_L_vals[avg_L_change==0]=np.nan

avg_bits_IH=np.nanmean(avg_I_vals,axis=0)
avg_bits_Lag=np.nanmean(avg_L_vals,axis=0)

LinRegSlope_H[LinRegSlope_H==0]=np.nan
LinRegSlope_L[LinRegSlope_L==0]=np.nan

avg_trend_IH=np.nanmean(LinRegSlope_H,axis=0)
avg_trend_L=np.nanmean(LinRegSlope_L,axis=0)
            