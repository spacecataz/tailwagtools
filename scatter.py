#import things
import matplotlib.pyplot as plt; plt.ion() #interactive matplotlib
import numpy as np
from tailwagtools import tailwag
import datetime as dt

#create an epoch
# epoch = dt.datetime(2001,8,19,20,0,0)
### Read list of events:
f = open('Event_Points.txt','r')
lines = f.readlines()
f.close()
bins = len(lines)
newlines = [None] * bins
data = [None] * bins
xdays = [None] * (bins-1)
days = [None] * (bins-1)
epoch = [None] * (bins-1)
t = [None] * (bins-1)
b89 = [None] * (bins-1)
b96 = [None] * (bins-1)
b01 = [None] * (bins-1)
loc = [None] * (bins-1)
for i in range(bins):
    newlines[i] = lines[i].replace('\n','')
    data[i] = newlines[i].split('\t')
for i in range(bins-1):
    xdays[i] = data[i+1][0]
    days[i] = xdays[i].replace('.000','')
    epoch[i] = dt.datetime.strptime(days[i], '%Y-%m-%d %H:%M:%S')
# day1 = data[1][0]

#load cluster data at desired time
for i in range(len(epoch)):
    data[i]  = tailwag.fetch_cluster_data(epoch[i]) # Use tspan here.
    t[i],b89[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T89')
    tdelta  = dt.timedelta(hours=6) #narrow down to +/- 6 hours
    loc[i]    = (t[i]>=epoch[i]-tdelta)&(t[i]<=epoch[i]+tdelta)
#get T89 data. see commentary on the errors i encountered below**
#Tricky indexing things! We like this!

bx_cluster = [None] * (bins-1)
t_cevent = [None] * (bins-1)
loc2 = [None] * (bins-1)
t_cluster = [None] * (bins-1)
#figure for 1 event
# fig, axs = plt.subplots(2)

#load and plot cluster data Bx at desired time
for i in range(len(epoch)):
    bx_cluster[i] = data[i]['b'][loc[i],0]
    t_cevent[i]   = data[i]['fgm_time'][loc[i]]
    loc2[i] = np.abs(bx_cluster[i])==np.abs(bx_cluster[i]).min() #crossing is min of |Bx|
    t_cluster[i] = t_cevent[i][loc2[i]]
#axs[0].plot(t_cevent,bx_cluster)
#pinpoint the exact crossing!
 #time of crossing
# axs[0].plot(t_clustebxr,[0],'b*', ms=5)

bx_t89 = [None] * (bins-1)
t_89event = [None] * (bins-1)
loc3 = [None] * (bins-1)
t_T89 = [None] * (bins-1)

#plot T89
for i in range(len(epoch)):
    bx_t89[i]  = b89[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
    t_89event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_89event,bx_t89)
#pinpoint the exact crossing!
    loc3[i]      = np.abs(bx_t89[i])==np.abs(bx_t89[i]).min() #crossing is min of |Bx|
    t_T89[i] = t_89event[i][loc3[i]] #time of crossing
# axs[0].plot(t_T89,[0],'b*', ms=5)

b96 = [None] * (bins-1)
bx_t96 = [None] * (bins-1)
t_96event = [None] * (bins-1)
loc4 = [None] * (bins-1)
t_T96 = [None] * (bins-1)

#get T96
for i in range(len(epoch)):
    t[i],b96[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T96')
    bx_t96[i]  = b96[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
    t_96event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_96event,bx_t96)
#pinpoint the exact crossing!
    loc4[i]      = np.abs(bx_t96[i])==np.abs(bx_t96[i]).min() #crossing is min of |Bx|
    t_T96[i] = t_96event[i][loc4[i]] #time of crossing
# axs[0].plot(t_T96,[0],'b*', ms=5)

b01 = [None] * (bins-1)
bx_t01 = [None] * (bins-1)
t_01event = [None] * (bins-1)
loc5 = [None] * (bins-1)
t_T01 = [None] * (bins-1)

#get T01STORM
for i in range(len(epoch)):
    t[i],b01[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T01STORM')
    bx_t01[i]  = b01[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
    t_01event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_01event,bx_t01)
#pinpoint the exact crossing!
    loc5[i]      = np.abs(bx_t01[i])==np.abs(bx_t01[i]).min() #crossing is min of |Bx|
    t_T01[i] = t_01event[i][loc5[i]] #time of crossing
# axs[0].plot(t_T01,[0],'b*', ms=5)

timediff_T89 = [None] * (bins-1)
for i in range(len(epoch)):
    timediff_T89[i] = abs(t_cluster[i]-t_T89[i])

t_cis = [None] * (bins-1)
loc6 = [None] * (bins-1)
index_cluster = [None] * (bins-1)
a = [None] * (bins-1)
b = [None] * (bins-1)
cluster_density = [None]

#plot densities
for i in range(len(epoch)):
    t_cis[i] = data[i]['cis_time']
    loc6[i]    = (t_cis[i]>=epoch[i]-tdelta)&(t_cis[i]<=epoch[i]+tdelta)
# axs[1].plot(t_cis[loc6],data['dens_h'][loc6])
# axs[1].plot(t_cis[loc6],data['dens_o'][loc6])

#*******************************

#identify plasma density at cluster time
for i in range(len(epoch)):
    index_cluster[i] = np.where(t_cis[i]<t_cluster[i])[0][-1]
    a[i] = t_cluster[i]-t_cis[i][index_cluster]
    b[i] = t_cis[i][index_cluster+1]-t_cluster[i]
    if a[i] > b[i]:
        index_cluster[i] = index_cluster[i]+1
    cluster_density[i] = data[i]['dens_h'][index_cluster] + data[i]['dens_o'][index_cluster]
# axs[1].plot(t_cis[index_cluster],data['dens_h'][index_cluster],'b*', ms=5)
# axs[1].plot(t_cis[index_cluster],data['dens_o'][index_cluster],'b*', ms=5)

#identify plasma density at T89 time
index_T89 = np.where(t_cis<t_T89)[0][-1]
a = t_T89-t_cis[index_T89]
b = t_cis[index_T89+1]-t_T89
if a > b:
    T89_cluster = index_T89+1
T89_density = data['dens_h'][index_T89] + data['dens_o'][index_T89]
# axs[1].plot(t_cis[index_T89],data['dens_h'][index_T89],'b*', ms=5)
# axs[1].plot(t_cis[index_T89],data['dens_o'][index_T89],'b*', ms=5)

#identify plasma density at T96 time
index_T96 = np.where(t_cis<t_T96)[0][-1]
a = t_T96-t_cis[index_T96]
b = t_cis[index_T96+1]-t_T96
if a > b:
    index_T96 = index_T96+1
T96_density = data['dens_h'][index_T96] + data['dens_o'][index_T96]
# axs[1].plot(t_cis[index_T96],data['dens_h'][index_T96],'b*', ms=5)
# axs[1].plot(t_cis[index_T96],data['dens_o'][index_T96],'b*', ms=5)

#identify plasma density at T01STORM time
index_T01STORM = np.where(t_cis<t_T01STORM)[0][-1]
a = t_T01STORM-t_cis[index_T01STORM]
b = t_cis[index_T01STORM+1]-t_T01STORM
if a > b:
    index_T01STORM = index_T01STORM+1
T01STORM_density = data['dens_h'][index_T01STORM] + data['dens_o'][index_T01STORM]
#axs[1].plot(t_cis[index_T01STORM],data['dens_h'][index_T01STORM],'b*', ms=5)
#axs[1].plot(t_cis[index_T01STORM],data['dens_o'][index_T01STORM],'b*', ms=5)
