# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 09:54:06 2019

@author: keckj
"""

#plot of max hourly precip rate versus mean daily rate for four NOAA gauges 
#AND one King county gauge

import matplotlib.pyplot as plt
import numpy as np
import os as os
import pandas as pd
os.chdir('C:/Users/keckj/Documents/GitHub/code/modules/')
import timeseriesstatistics as TSS
import scipy.interpolate as interpolate
import matplotlib.dates as mdates
import datetime
os.chdir('C:/Users/keckj/Documents/GitHub/code/modules/')
import HydroGeoFunctions as hgf
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)




#compare change in SWE at time of storm
#load SNOTEL data

os.chdir('D:/DNR/Culverts/tempdata')

yrs = np.linspace(1990,2013,(2013-1990+1))

#concat all years into a single data set
for i in yrs:
    #i=1990
    tdat = pd.read_csv('817_STAND_WATERYEAR='+str(int(i))+'.csv', sep=',', header=1)
    if i == yrs[0]:
        DAT = tdat
    else:
        DAT = pd.concat([DAT,tdat],axis=0)



DAT['WTEQdf'] = DAT['WTEQ.I-1 (in) '].diff()
DAT['DTdate'] = pd.to_datetime(DAT['Date'])
DAT = DAT.set_index('DTdate')

os.chdir('D:/UW_PhD/PreeventsProject/Hydrometeorology and Floods/HydrologyObservations')

#Sauk mountain gauges
#data = pd.read_csv('SaukPrecip_HourlyObservations.csv')
#
##NOAA hourly pricp data: only hours of P recorded, all missing data is zero
##'STEHEKIN 4 NW WA US' misssing data, using other three gages
#nms = ['DARRINGTON RANGER STATION WA US','DIABLO DAM WA US','MARBLEMOUNT RANGER STATION WA US','BURLINGTON WA US']
#Sdat = {}
#for i in nms:
#    print(i[0:4])
#    
#    dat = data.loc[data['STATION_NAME']==i].copy()
#    dat['HPCP'].loc[dat['HPCP'] == 25399.75] = 0 #np.nan
#    dat['DATE'] = pd.to_datetime(dat['DATE']) #convert date column to datetime64 fromat
#    dat = dat.set_index('DATE', drop=False) # set datetime64 as index
#    dat = dat.resample('H').asfreq() #resample data set as hourly so that all hours without data are nan
#    dat['HPCP'] = dat['HPCP'].fillna(0) #fill nans with 0
#    dat['HPCP'].loc[dat['HPCP'] > 30] = 0 #fill nans with 0
#    dat['DPCP'] = dat['HPCP'].resample('D').mean() # create mean daily column
#    dat['DPCP']= dat['DPCP'].fillna(method='ffill')
#    Sdat[i[0:4]] = dat
#
#
##date begin 1984-10-01. date end: 2012-09-30 for DIAB, MARB, 2008-09-30 for DARR 
#Snms = ['DARR', 'DIAB', 'MARB','BURL']
#dbl = ['1984-10-01','1984-10-01','1984-10-01']
#denl = ['2008-09-30','2012-09-30','2012-09-30'] 
#
#sta = 0
#db = '1990'#dbl[sta]
#de = '2013-09-30'#denl[sta]
#dat = Sdat[Snms[sta]]

'''
###king county data
'''
db = pd.to_datetime('2000')#dbl[sta]
de = pd.to_datetime('2013-09-30')#denl[sta]
dat = pd.read_csv('Hydrology_DCJKI_KingCountySEREStation.csv',index_col=False)


dat['Date'] = pd.to_datetime(dat['Collect Date (local)']) #convert date column to datetime64 fromat
dat = dat.set_index('Date', drop=False) # set datetime64 as index
#dat = dat.resample('H').asfreq() #resample data set as hourly so that all hours without data are nan
dat['Precipitation (inches)'] = dat['Precipitation (inches)'].fillna(0) #fill nans with 0
dat['HPCP'] = dat['Precipitation (inches)']*25.4 #np.nan
dat['DPCP'] = dat['HPCP'].resample('D').mean() # create mean daily column
dat['DPCP']= dat['DPCP'].fillna(method='ffill')

yb = dat.index.year.min()
ye = 2012#dat.index.year.max()
yn = ye-yb+1
yrs = np.linspace(yb,ye,yn)

#visually examine each year of data
for i in yrs:
    yrp = int(i)
    sdb= pd.to_datetime(str(yrp)+'-10-01')
    sde = pd.to_datetime(str(yrp+1)+ '-09-01')
    cc = dat[sdb:sde]['HPCP']
    dd = dat[sdb:sde]['HPCP'].cumsum()
    fig,ax = plt.subplots(figsize=(18,12))
    plt.plot(dd.index,dd.values)
    plt.xlim([sdb,sde])
    plt.ylim([dd[sdb:sde].min(),dd[sdb:sde].max()])
    plt.show()
    
    fig,ax = plt.subplots(figsize=(18,3))
    plt.plot(dat.index,dat['HPCP'])
    plt.xlim([db,de])
    plt.ylim([dat['HPCP'][sdb:sde].min()*.96,dat['HPCP'][sdb:sde].max()*1.02])
    plt.show()




fd1 = dat['HPCP']
fd2 = dat['DPCP']
fd1_am = fd1.resample('Y').max()
fd2_am = fd2.resample('Y').max()
fd1_t = fd1[db:de] #already truncated
fd2_t = fd2[db:de] #already truncated


# water year annual max record
df = pd.DataFrame(fd1_t)
WYs_basin = TSS.WaterYear_Metric(df,metric = 'max',qval=None)
WY_sta = WYs_basin[0]
WY_sta_t  = WY_sta[db:de]


#partial duration series of max
bS = {'Inst':[fd1]}
dnm = ['Inst']
out = hgf.PartialDurationSeries(dnm,bS,Stype='PDS',RI=[5],events = False) #use the entire record to compute RI

pds = out[0]['Inst']


pds_t = pds[str(db.year)+'-10-02':de.strftime('%Y-%m-%d')] # extract all events in pds that can be paired with the SWE record
WY_sta_t = pds_t.reset_index()
WY_sta_t.columns = ['Date','P [mm/hr]','RI']

WTEQdf = []
#determine change SWE for day of mean peak
for index, row in WY_sta_t.iterrows(): 
    aa = row['Date'].strftime('%Y-%m-%d')
    WTEQdf.append(DAT[DAT['Date'] == aa]['WTEQdf'].values[0])

WTEQdfd = pd.DataFrame(WTEQdf,index = WY_sta_t.index, columns = ['SWE_dif [in]'])
WY_sta_t = pd.concat([WY_sta_t,WTEQdfd],axis=1)
WY_sta_t = WY_sta_t.set_index('Date',drop =False)



#beacuase partial duration series do not MATCH, MATCH daily peak to hourly to create WY_sta_t2

date_c = []
dm_c = []
for index, row in WY_sta_t.iterrows():
    print(row['Date'].day)
    ad = row['Date'].day
    aa = row['Date'].strftime('%Y-%m-%d')
    
    #day before, of and after date
    dt1 = pd.to_datetime(aa)+ datetime.timedelta(days=-1)
    dt2 = pd.to_datetime(aa)+ datetime.timedelta(days=0)
    dt3 = pd.to_datetime(aa)+ datetime.timedelta(days=1)
    dtd = pd.DataFrame([dt1,dt2,dt3])

    
    #day before, of and after date value
    D1 = fd2[dt1]
    D2 = fd2[dt2]
    D3 = fd2[dt3]
    dd = np.array([D1,D2,D3])
    
    maxdm = dd[dd == dd.max()][0]
    maxdmdt = dtd[0][dd == dd.max()].iloc[0] 
    dm_c.append(maxdm)
    date_c.append(maxdmdt)
    
WY_sta_t2 = pd.DataFrame(list(zip(dm_c,date_c)),columns = ['P [mm/hr]','Date_m'])
WY_sta_t2.index = WY_sta_t.index



#determine change WY for day of mean peak
WTEQdf2 = []
for index, row in WY_sta_t2.iterrows(): 
    aa = row['Date_m'].strftime('%Y-%m-%d')
    WTEQdf2.append(DAT[DAT['Date'] == aa]['WTEQdf'].values[0])

WTEQdfd2 = pd.DataFrame(WTEQdf2,index = WY_sta_t2.index, columns = ['SWE_dif [in]'])

WY_sta_t2 = pd.concat([WY_sta_t2,WTEQdfd2],axis=1)


#add daily mean flow value and ratio of instantaneous to daily mean to WY_sta_t

WY_sta_t['Dmean P [mm/hr]'] = WY_sta_t2['P [mm/hr]']
WY_sta_t['flow ratio'] = WY_sta_t['P [mm/hr]']/WY_sta_t2['P [mm/hr]']




fig, ax = plt.subplots(figsize=(8, 8))
c1 = True
c2= True
cc = 0
for index, row in WY_sta_t2.iterrows():
    x = row['P [mm/hr]']
    y = WY_sta_t['P [mm/hr]'].iloc[cc]
    if (row['SWE_dif [in]'] > 0) & (y>x): 
        if c1:
            plt.plot(x,y,'bo', label = 'Snow accumulation')
            c1 = False
        else:
            plt.plot(x,y,'bo', label = '_nolabel_')
    elif (row['SWE_dif [in]'] <= 0) & (y>x): 
        if c2:
            plt.plot(x,y,'ro', label = 'Snow melt/no change')
            c2 = False
        else:
            plt.plot(x,y,'ro', label = '_nolabel_')
            
    cc +=1
    
plt.plot([0,10e7],[0,10e7],'k-')
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.xticks(rotation=30)
plt.xlabel('max mean daily P [mm/hr]',fontsize = 20)
plt.ylabel('max mean hourly P [mm/hr]',fontsize = 20)
plt.xlim([0,10])
plt.ylim([0,max(fd1_am.values)*1.1])
plt.legend(fontsize = 12, loc = 'lower right')
#plt.savefig('ObservedPHourlyvsDailyMeanRate_'+Snms[sta]+'.jpg', dpi = 300, bbox_inches='tight')
plt.savefig('ObservedPHourlyvsDailyMeanRate_KingCounty.jpg', dpi = 300, bbox_inches='tight')
plt.show()

