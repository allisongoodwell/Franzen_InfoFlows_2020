# -*- coding: utf-8 -*-
"""
Updated 6/5/2020

This is modified from Sam Franzen's original code, to input flow and precip data
compute information measures for entire time window
primary metrics: mutual information, dominant lag time, threhsold, specific information values

Inputs: CPC gage-based gridded data for CO Headwaters HUC4 basin region, USGS flow rate data

#functions: comput_info_measures.py and compute_icrit.py for IT computations

Outputs: 

@author: Sam, Allison, Mozhgan
"""

import pickle
import numpy as np
import pandas as pd
import time
from itertools import chain
from compute_icrit import compute_icrit as ci
from compute_info_measures import compute_info_measures as cim

start = time.time()

#designate bins for P, xbin, and Q, ybin (change it from 5 to 17)  
xbin = 2
ybin = 5

#precipitation thresholds (mm per day)
thresholds = [0.3,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

#designate max lag
delay = 7

start_year1 = 1953
end_year1 = 2016

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


#load streamflow data (It could be "COHW flow rates" & "difference between COHW and Gunnison flow rates")
dfsf1 = pd.read_csv('COHW Gage Data 1952-2017.csv', usecols=[4])
dfsf1.rename(columns={'difference': 'flow'}, inplace=True)
dfsf1_list = dfsf1['flow']
dfsf1_list = dfsf1_list.tolist()
dfsf1_series = pd.Series(dfsf1_list)
dfsf1_array = np.array(dfsf1_series)

#load gridded rainfall data
f_myfile = open('CO_rainfall_data.pickle','rb')
Co_pptdata = pickle.load(f_myfile) #ppt_long_list = rainfall data [0] and percentiles [1]
lat = pickle.load(f_myfile)
lon = pickle.load(f_myfile)
f_myfile.close()

dfppt = pd.DataFrame(Co_pptdata)


yearvect = range(start_year1,end_year1)

#create lists and arrays for data storage
winterplist = []
springplist = []
summerplist = []
fallplist = []

winterflist = []
springflist = []
summerflist = []
fallflist = []

winterlen = np.sum([31,31,28]) #DJF
springlen = np.sum([31,30,31]) #MAM
summerlen = np.sum([30, 31, 31]) #JJA
falllen = np.sum([30, 31, 30]) #SON

#breakpoints between years within a 5-year period
winter_breakpoints = [winterlen*i for i in range(1,3)]
spring_breakpoints = [springlen*i for i in range(1,3)]
summer_breakpoints = [summerlen*i for i in range(1,3)]
fall_breakpoints = [falllen*i for i in range(1,3)]

winterdata = np.zeros((3,len(Co_pptdata)))
springdata = np.zeros((3,len(Co_pptdata)))
summerdata = np.zeros((3,len(Co_pptdata)))
falldata = np.zeros((3,len(Co_pptdata)))

winterlags = np.zeros((len(Co_pptdata)))
springlags = np.zeros((len(Co_pptdata)))
summerlags = np.zeros((len(Co_pptdata)))
falllags = np.zeros((len(Co_pptdata)))

winterIp = np.zeros((3,len(Co_pptdata),len(thresholds),delay+1))
springIp = np.zeros((3,len(Co_pptdata),len(thresholds),delay+1))
summerIp = np.zeros((3,len(Co_pptdata),len(thresholds),delay+1))
fallIp = np.zeros((3,len(Co_pptdata),len(thresholds),delay+1))

winterHQ = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
springHQ = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
summerHQ = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
fallHQ = np.zeros((len(Co_pptdata),len(thresholds),delay+1))

winterSI = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
springSI = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
summerSI = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
fallSI = np.zeros((len(Co_pptdata),len(thresholds),delay+1))

winterp1 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
springp1 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
summerp1 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
fallp1 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))

winterp0 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
springp0 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
summerp0 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))
fallp0 = np.zeros((len(Co_pptdata),len(thresholds),delay+1))

wintermaxI = np.zeros((len(Co_pptdata)))
springmaxI = np.zeros((len(Co_pptdata)))
summermaxI = np.zeros((len(Co_pptdata)))
fallmaxI = np.zeros((len(Co_pptdata)))

wintermaxIrain = np.zeros((len(Co_pptdata)))
springmaxIrain = np.zeros((len(Co_pptdata)))
summermaxIrain = np.zeros((len(Co_pptdata)))
fallmaxIrain = np.zeros((len(Co_pptdata)))

wintermaxIdry = np.zeros((len(Co_pptdata)))
springmaxIdry = np.zeros((len(Co_pptdata)))
summermaxIdry = np.zeros((len(Co_pptdata)))
fallmaxIdry = np.zeros((len(Co_pptdata)))

wintermaxSI = np.zeros((len(Co_pptdata)))
springmaxSI = np.zeros((len(Co_pptdata)))
summermaxSI = np.zeros((len(Co_pptdata)))
fallmaxSI = np.zeros((len(Co_pptdata)))

wintermaxp1 = np.zeros((len(Co_pptdata)))
springmaxp1 = np.zeros((len(Co_pptdata)))
summermaxp1 = np.zeros((len(Co_pptdata)))
fallmaxp1 = np.zeros((len(Co_pptdata)))

wintermaxp0 = np.zeros((len(Co_pptdata)))
springmaxp0 = np.zeros((len(Co_pptdata)))
summermaxp0 = np.zeros((len(Co_pptdata)))
fallmaxp0 = np.zeros((len(Co_pptdata)))

wintermaxHQ = np.zeros((len(Co_pptdata)))
springmaxHQ = np.zeros((len(Co_pptdata)))
summermaxHQ = np.zeros((len(Co_pptdata)))
fallmaxHQ = np.zeros((len(Co_pptdata)))

wintermaxthresh = np.zeros((len(Co_pptdata)))
springmaxthresh = np.zeros((len(Co_pptdata)))
summermaxthresh = np.zeros((len(Co_pptdata)))
fallmaxthresh = np.zeros((len(Co_pptdata)))

wintermaxL = np.zeros((len(Co_pptdata)))
springmaxL = np.zeros((len(Co_pptdata)))
summermaxL = np.zeros((len(Co_pptdata)))
fallmaxL = np.zeros((len(Co_pptdata)))

winterdata = np.zeros((3,len(Co_pptdata)))
springdata = np.zeros((3,len(Co_pptdata)))
summerdata = np.zeros((3,len(Co_pptdata)))
falldata = np.zeros((3,len(Co_pptdata)))

winterdatathresh = np.zeros((3,len(Co_pptdata)))
springdatathresh = np.zeros((3,len(Co_pptdata)))
summerdatathresh = np.zeros((3,len(Co_pptdata)))
falldatathresh = np.zeros((3,len(Co_pptdata)))

winterdataL = np.zeros((3,len(Co_pptdata)))
springdataL = np.zeros((3,len(Co_pptdata)))
summerdataL = np.zeros((3,len(Co_pptdata)))
falldataL = np.zeros((3,len(Co_pptdata)))

#load ppt data, specific station
for pg in range (0, len(Co_pptdata)):
    pptr8 = []
    
    pptr8 = list(chain(*Co_pptdata[pg]))
    del pptr8[0:1461]
    
    dfpptr8 = pd.DataFrame(pptr8)

    dfppsf = pd.concat([dfpptr8, dfsf1], axis = 1)
            
    for thvect,th in enumerate(thresholds):
        
        plimit = th
        
        bi_ppt_list = []
    
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
        
        start_year = 1956
        end_year = 2017
        
    
        for y in range(start_year,end_year):

            for s in range(1,5):
              
                if s>1:  #spring, summer, fall get data for window around year
                    
                    df_short = dfcc[np.logical_and(dfcc['season']==s,dfcc['year']==y)]
                    plist = df_short['ppt8'].tolist()
                    flist = df_short['flow'].tolist()
    
                    if s == 2:
                        springplist.extend(plist)
                        springflist.extend(flist)
        
                    if s == 3:
                        summerplist.extend(plist)
                        summerflist.extend(flist)
                        
                    if s == 4:
                        fallplist.extend(plist)
                        fallflist.extend(flist)
                    
                else: #winter is awkward (want Dec of past year and Jan-Feb of this year) get data for window around year
                                               
                    df_short = dfcc[np.logical_or(np.logical_and(dfcc['month']<3,dfcc['year']==y),np.logical_and(dfcc['month']==12,dfcc['year']==y-1))]
                    plist = df_short['ppt8'].tolist()
                    flist = df_short['flow'].tolist()
                    
                    winterplist.extend(plist)
                    winterflist.extend(flist)
                
        
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
                
            
            for d in range(0, delay+1):
                
                pptlist_new = pptlist[0:len(pptlist)-d] 
                sflist_new = sflist[d:len(sflist)]
                
                #AEG: here is the strange part where we have to leave out certain points...(breakpoints)
                if d>0:
                    for b in sorted(breakpoints,reverse=True):
                        pptlist_new[(b-d+1):b]=[]
                        sflist_new[(b-d+1):b]=[]
                        
                
                #create pdf of data and calculate information
                hist_array = np.array((pptlist_new, sflist_new)).T
                pdf, edges = np.histogramdd(hist_array, (xbin, ybin))
                pdf = pdf/np.sum(pdf)
                I_dict = cim(pdf)
                Icrit_dict = ci(pptlist_new, sflist_new)
                
                #Check for significance and save data to vectors for each season
                if t == 1:
                    a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                    if a > 0:
                        winterIp[0,pg,thvect,d]=I_dict['Ix0']
                        winterIp[1,pg,thvect,d]=I_dict['Ix1']
                        winterIp[2,pg,thvect,d]=I_dict['Ix1x2']
                        winterSI[pg,thvect,d]=I_dict['SI']
                        winterp1[pg,thvect,d]=I_dict['m1']
                        winterp0[pg,thvect,d]=I_dict['m0']
                        
                        winterHQ[pg,thvect,d]=I_dict['Hx2']
                    
                
                if t == 2:
                    a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                    if a > 0:
                        springIp[0,pg,thvect,d]=I_dict['Ix0']
                        springIp[1,pg,thvect,d]=I_dict['Ix1']
                        springIp[2,pg,thvect,d]=I_dict['Ix1x2']
                        springSI[pg,thvect,d]=I_dict['SI']
                        springp1[pg,thvect,d]=I_dict['m1']
                        springp0[pg,thvect,d]=I_dict['m0']
                        
                        springHQ[pg,thvect,d]=I_dict['Hx2']
                
                if t == 3:
                    a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                    if a > 0:
                        summerIp[0,pg,thvect,d]=I_dict['Ix0']
                        summerIp[1,pg,thvect,d]=I_dict['Ix1']
                        summerIp[2,pg,thvect,d]=I_dict['Ix1x2']
                        summerSI[pg,thvect,d]=I_dict['SI']
                        summerp1[pg,thvect,d]=I_dict['m1']
                        summerp0[pg,thvect,d]=I_dict['m0']
                        
                        summerHQ[pg,thvect,d]=I_dict['Hx2']
                
                if t == 4:
                    a = I_dict['Ix1x2'] - Icrit_dict['Icrit']
                    if a > 0:
                        fallIp[0,pg,thvect,d]=I_dict['Ix0']
                        fallIp[1,pg,thvect,d]=I_dict['Ix1']
                        fallIp[2,pg,thvect,d]=I_dict['Ix1x2']
                        fallSI[pg,thvect,d]=I_dict['SI']
                        fallp1[pg,thvect,d]=I_dict['m1']
                        fallp0[pg,thvect,d]=I_dict['m0']
                        
                        fallHQ[pg,thvect,d]=I_dict['Hx2']
        
         
#%%
# find maximum mutual information and the threshold and lag associated with the maximum MI for each gage and season
for pgi in range (0,43):

    #winter
    reg = int(2)
    wintermaxI[pgi]= np.max(winterIp[reg,pgi].flatten(), axis=0)
    maxthresh = int(np.ceil(np.argmax(winterIp[reg,pgi])/8))
    if maxthresh == 0:
        maxthresh = 0
    else: 
        maxthresh = maxthresh-1
    if maxthresh == 0:
        maxthresh1 = 0.3
    else: 
        maxthresh1 = maxthresh
    wintermaxthresh[pgi]=maxthresh1
    maxlag = int(np.argmax(winterIp[reg,pgi,maxthresh]))
    wintermaxL[pgi] = maxlag
    reg = int(1)
    wintermaxIrain[pgi] = winterIp[reg,pgi,maxthresh,maxlag]
    reg = int(0)
    wintermaxIdry[pgi] = winterIp[reg,pgi,maxthresh,maxlag]
    wintermaxSI[pgi] = winterSI[pgi,maxthresh,maxlag]
    wintermaxp1[pgi] = winterp1[pgi,maxthresh,maxlag]
    wintermaxp0[pgi] = winterp0[pgi,maxthresh,maxlag]
    wintermaxHQ[pgi] = winterHQ[pgi,maxthresh,maxlag]
    singleyearlist = list(range(start_year1,end_year1+1,3))
    x = np.asarray(singleyearlist)
    y = wintermaxI[pgi]
        
    #spring
    reg = 2
    springmaxI[pgi]= np.max(springIp[reg,pgi].flatten(), axis=0)
    maxthresh=int(np.ceil(np.argmax(springIp[reg,pgi])/8))
    if maxthresh == 0:
        maxthresh = 0
    else: 
        maxthresh = maxthresh-1
    if maxthresh == 0:
        maxthresh1 = 0.3
    else: 
        maxthresh1 = maxthresh
    springmaxthresh[pgi]=maxthresh1
    maxlag=int(np.argmax(springIp[reg,pgi,maxthresh]))
    springmaxL[pgi]=maxlag
    reg = 1
    springmaxIrain[pgi] = springIp[reg,pgi,maxthresh,maxlag]
    reg = 0
    springmaxIdry[pgi] = springIp[reg,pgi,maxthresh,maxlag]
    springmaxSI[pgi] = springSI[pgi,maxthresh,maxlag]
    springmaxp1[pgi] = springp1[pgi,maxthresh,maxlag]
    springmaxp0[pgi] = springp0[pgi,maxthresh,maxlag]
    springmaxHQ[pgi] = springHQ[pgi,maxthresh,maxlag]
    singleyearlist = list(range(start_year1,end_year1+1,3))
    x = np.asarray(singleyearlist)
    y = springmaxI[pgi]
    
    #summer
    reg = 2
    summermaxI[pgi]= np.max(summerIp[reg,pgi].flatten(), axis=0)
    maxthresh=int(np.ceil(np.argmax(summerIp[reg,pgi])/8))
    if maxthresh == 0:
        maxthresh = 0
    else: 
        maxthresh = maxthresh-1
    if maxthresh == 0:
        maxthresh1 = 0.3
    else: 
        maxthresh1 = maxthresh
    summermaxthresh[pgi]=maxthresh1
    maxlag=int(np.argmax(summerIp[reg,pgi,maxthresh]))
    summermaxL[pgi]=maxlag
    reg = 1
    summermaxIrain[pgi] = summerIp[reg,pgi,maxthresh,maxlag]
    reg = 0
    summermaxIdry[pgi] = summerIp[reg,pgi,maxthresh,maxlag]
    summermaxSI[pgi] = summerSI[pgi,maxthresh,maxlag]
    summermaxp1[pgi] = summerp1[pgi,maxthresh,maxlag]
    summermaxp0[pgi] = summerp0[pgi,maxthresh,maxlag]
    summermaxHQ[pgi] = summerHQ[pgi,maxthresh,maxlag]
    singleyearlist = list(range(start_year1,end_year1+1,3))
    x = np.asarray(singleyearlist)
    y = summermaxI[pgi]
    
    #fall
    reg=2
    fallmaxI[pgi]= np.max(fallIp[reg,pgi].flatten(), axis=0)
    maxthreshyr=int(np.ceil(np.argmax(fallIp[reg,pgi])/8))
    if maxthresh == 0:
        maxthresh = 0
    else: 
        maxthresh = maxthresh-1
    if maxthresh == 0:
        maxthresh1 = 0.3
    else: 
        maxthresh1 = maxthresh
    fallmaxthresh[pgi]=maxthresh1
    maxlag=int(np.argmax(fallIp[reg,pgi,maxthresh]))
    fallmaxL[pgi]=maxlag
    reg = 1
    fallmaxIrain[pgi] = fallIp[reg,pgi,maxthresh,maxlag]
    reg = 0
    fallmaxIdry[pgi] = fallIp[reg,pgi,maxthresh,maxlag]
    fallmaxSI[pgi] = fallSI[pgi,maxthresh,maxlag]
    fallmaxp1[pgi] = fallp1[pgi,maxthresh,maxlag]
    fallmaxp0[pgi] = fallp0[pgi,maxthresh,maxlag]
    fallmaxHQ[pgi] = fallHQ[pgi,maxthresh,maxlag]
    singleyearlist = list(range(start_year1,end_year1+1,3))
    x = np.asarray(singleyearlist)
    y = fallmaxI[pgi]
  
#%%
#save results for overall analysis in an excel file
lon1 = [-x for x in lon]

dflat=pd.DataFrame(lat, columns = ['lat'])
dflon1=pd.DataFrame(lon1, columns = ['lon1'])

dfwintermaxI = pd.DataFrame(wintermaxI, columns = ['winterI'])  
dfspringmaxI = pd.DataFrame(springmaxI, columns = ['springI'])
dfsummermaxI = pd.DataFrame(summermaxI, columns = ['summerI'])
dffallmaxI = pd.DataFrame(fallmaxI, columns = ['fallI'])

dfwintermaxIrain = pd.DataFrame(wintermaxIrain, columns = ['winterIrain'])  
dfspringmaxIrain = pd.DataFrame(springmaxIrain, columns = ['springIrain'])
dfsummermaxIrain = pd.DataFrame(summermaxIrain, columns = ['summerIrain'])
dffallmaxIrain = pd.DataFrame(fallmaxIrain, columns = ['fallIrain'])

dfwintermaxIdry = pd.DataFrame(wintermaxIdry, columns = ['winterIdry'])  
dfspringmaxIdry = pd.DataFrame(springmaxIdry, columns = ['springIdry'])
dfsummermaxIdry = pd.DataFrame(summermaxIdry, columns = ['summerIdry'])
dffallmaxIdry = pd.DataFrame(fallmaxIdry, columns = ['fallIdry'])
            
dfwinterthresh = pd.DataFrame(wintermaxthresh, columns = ['winterT'])  
dfspringthresh = pd.DataFrame(springmaxthresh, columns = ['springT'])
dfsummerthresh = pd.DataFrame(summermaxthresh, columns = ['summerT'])
dffallthresh = pd.DataFrame(fallmaxthresh, columns = ['fallT'])

dfwinterL = pd.DataFrame(wintermaxL, columns = ['winterL'])  
dfspringL = pd.DataFrame(springmaxL, columns = ['springL'])
dfsummerL = pd.DataFrame(summermaxL, columns = ['summerL'])
dffallL = pd.DataFrame(fallmaxL, columns = ['fallL'])

dfwinterSI = pd.DataFrame(wintermaxSI, columns = ['winterSI'])  
dfspringSI = pd.DataFrame(springmaxSI, columns = ['springSI'])
dfsummerSI = pd.DataFrame(summermaxSI, columns = ['summerSI'])
dffallSI = pd.DataFrame(fallmaxSI, columns = ['fallSI'])

dfwintermaxp0 = pd.DataFrame(wintermaxp0, columns = ['winterp0'])  
dfspringmaxp0 = pd.DataFrame(springmaxp0, columns = ['springp0'])
dfsummermaxp0 = pd.DataFrame(summermaxp0, columns = ['summerp0'])
dffallmaxp0 = pd.DataFrame(fallmaxp0, columns = ['fallp0'])

dfwintermaxp1 = pd.DataFrame(wintermaxp1, columns = ['winterp1'])  
dfspringmaxp1 = pd.DataFrame(springmaxp1, columns = ['springp1'])
dfsummermaxp1 = pd.DataFrame(summermaxp1, columns = ['summerp1'])
dffallmaxp1 = pd.DataFrame(fallmaxp1, columns = ['fallp1'])

dfwintermaxHQ = pd.DataFrame(wintermaxHQ, columns = ['winterHQ'])  
dfspringmaxHQ = pd.DataFrame(springmaxHQ, columns = ['springHQ'])
dfsummermaxHQ = pd.DataFrame(summermaxHQ, columns = ['summerHQ'])
dffallmaxHQ = pd.DataFrame(fallmaxHQ, columns = ['fallHQ'])

df60results = pd.concat([dflat, dflon1, dfwintermaxI, dfspringmaxI, dfsummermaxI, dffallmaxI, dfwintermaxHQ, dfspringmaxHQ, dfsummermaxHQ, dffallmaxHQ, dfwintermaxIrain, dfspringmaxIrain, dfsummermaxIrain, dffallmaxIrain, \
                       dfwintermaxIdry, dfspringmaxIdry, dfsummermaxIdry, dffallmaxIdry, dfwinterthresh, dfspringthresh, dfsummerthresh, dffallthresh, \
                       dfwinterL, dfspringL, dfsummerL, dffallL, dfwinterSI, dfspringSI, dfsummerSI, dffallSI, \
                       dfwintermaxp0, dfspringmaxp0, dfsummermaxp0, dffallmaxp0, dfwintermaxp1, dfspringmaxp1, dfsummermaxp1, dffallmaxp1], axis = 1)

df60results.to_csv('Overall_Results.csv')

end = time.time()
print ((end-start)/60)   




