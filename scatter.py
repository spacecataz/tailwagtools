#import things
import matplotlib.pyplot as plt; plt.ion() #interactive matplotlib
import numpy as np
from tailwagtools import tailwag
import datetime as dt
import scipy
from scipy import stats
import math

#create an epoch
# epoch = dt.datetime(2001,8,19,20,0,0)
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
    data[i]  = tailwag.fetch_cluster_data(epoch[i])
    t[i],b89[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T89')
    tdelta  = dt.timedelta(hours=6) #narrow down to +/- 6 hours
    loc[i]    = (t[i]>=epoch[i]-tdelta)&(t[i]<=epoch[i]+tdelta)


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
timediff_T89 = [None] * (bins-1)


#plot T89
for i in range(len(epoch)):
    bx_t89[i]  = b89[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
    t_89event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_89event,bx_t89)
#pinpoint the exact crossing!
    loc3[i]      = np.abs(bx_t89[i])==np.abs(bx_t89[i]).min() #crossing is min of |Bx|
    t_T89[i] = t_89event[i][loc3[i]] #time of crossing
    timediff_T89[i] = abs(t_cluster[i]-t_T89[i])
# axs[0].plot(t_T89,[0],'b*', ms=5)

# --------------------------------------
b96 = [None] * (bins-1)
bx_t96 = [None] * (bins-1)
t_96event = [None] * (bins-1)
loc4 = [None] * (bins-1)
t_T96 = [None] * (bins-1)
timediff_T96 = [None] * (bins-1)

#get T96
for i in range(len(epoch)):
    t[i],b96[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T96')
    bx_t96[i]  = b96[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
    t_96event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_96event,bx_t96)
#pinpoint the exact crossing!
    loc4[i]      = np.abs(bx_t96[i])==np.abs(bx_t96[i]).min() #crossing is min of |Bx|
    t_T96[i] = t_96event[i][loc4[i]] #time of crossing
    timediff_T96[i] = abs(t_cluster[i]-t_T96[i])
# axs[0].plot(t_T96,[0],'b*', ms=5)

b01 = [None] * (bins-1)
bx_t01 = [None] * (bins-1)
t_01event = [None] * (bins-1)
loc5 = [None] * (bins-1)
t_T01 = [None] * (bins-1)
timediff_T01 = [None] * (bins-1)


#get T01STORM
for i in range(len(epoch)):
    try:
        t[i],b01[i] = tailwag.gen_sat_tsyg(data[i],extMag = 'T01STORM')
        bx_t01[i]  = b01[i][loc[i],0] #reduce B to only Bx during interested time, shape (10780,)
        t_01event[i] = t[i][loc[i]] #define time of event
# axs[0].plot(t_01event,bx_t01)
#pinpoint the exact crossing!
        loc5[i]      = np.abs(bx_t01[i])==np.abs(bx_t01[i]).min() #crossing is min of |Bx|
        t_T01[i] = t_01event[i][loc5[i]] #time of crossing
        timediff_T01[i] = abs(t_cluster[i]-t_T01[i])
    #except TypeError: WORK ON THIS!
        #t_T01[i] = #change this shape somehow, ValueError: operands could not be broadcast together with shapes (2,) (0,)

# axs[0].plot(t_T01,[0],'b*', ms=5)
*************************************************
# plasma calculations!

t_cis = [None] * (bins-1)
t_plasma = [None] * (bins-1)
loc6 = [None] * (bins-1)
index_cluster = [None] * (bins-1)
a = [None] * (bins-1)
b = [None] * (bins-1)
cluster_density = [None]
h_dens = [None] * (bins-1)
o_dens = [None] * (bins-1)

#plot densities
for i in range(len(epoch)):
    t_cis[i] = data[i]['cis_time']
    #loc6[i]    = (t_cis[i]>=epoch[i]-tdelta)&(t_cis[i]<=epoch[i]+tdelta)
    #t_plasma[i] = t_cis[i][loc6[i]]
    h_dens[i] = data[i]['dens_h']
    o_dens[i] = data[i]['dens_o']
# axs[1].plot(t_cis[loc6],data['dens_h'][loc6])
# axs[1].plot(t_cis[loc6],data['dens_o'][loc6])

# WEIRD EVENT: t_cluster[63].shape = (2,)

for i in range(len(epoch)):
    try:
        index_cluster[i] = np.where(t_cis[i]<t_cluster[i])[0][-1]
        a[i] = t_cluster[i]-t_cis[i][index_cluster[i]]
        b[i] = t_cis[i][index_cluster[i]+1]-t_cluster[i]
        if a[i] > b[i]:
            index_cluster[i] = index_cluster[i]+1
    except IndexError:
        index_cluster[i] = [None]


loc7 = [None] * (bins-1)
cluster_after = [None] * (bins-1)
loc8 = [None] * (bins-1)
cluster_before = [None] * (bins-1)


for i in range(len(epoch)):
    try:
        loc7[i] = t_cis[i]>t_cis[i][index_cluster[i]]
        cluster_after[i] = h_dens[i][loc7[i]].mean() + o_dens[i][loc7[i]].mean()
    except IndexError:
        print(np.array(loc7[i]).shape)
        print(loc7[i])
        cluster_after[i] = [None]
        print(cluster_after[i])
for i in range(len(epoch)):
    try:
        loc8[i] = t_cis[i]<t_cis[i][index_cluster[i]]
        cluster_before[i] = h_dens[i][loc8[i]].mean() + o_dens[i][loc8[i]].mean()
    except IndexError:
        print(np.array(loc8[i]).shape)
        print(loc8[i])
        cluster_before[i] = [None]
        print(cluster_before[i])

cluster_diff = [None] * (bins-1)

for i in range(len(epoch)):
    try:
        cluster_diff[i] = abs(cluster_after[i]-cluster_before[i])
    except TypeError:
        print(cluster_after[i])
        print(cluster_before[i])
        cluster_diff[i] = [None]
        

# axs[1].plot(t_cis[index_cluster],data['dens_h'][index_cluster],'b*', ms=5)
# axs[1].plot(t_cis[index_cluster],data['dens_o'][index_cluster],'b*', ms=5)

#identify plasma density at T89 time
index_T89 = [None] * (bins-1)
a = [None] * (bins-1)
b = [None] * (bins-1)

for i in range(len(epoch)):
    try:
        index_T89[i] = np.where(t_cis[i]<t_T89[i])[0][-1]
        a[i] = t_T89[i]-t_cis[i][index_T89[i]]
        b[i] = t_cis[i][index_T89[i]+1]-t_T89[i]
        if a[i] > b[i]:
            index_T89[i] = index_T89[i]+1
    except IndexError:
        index_T89[i] = [None]

loc9 = [None] * (bins-1)
T89_after = [None] * (bins-1)
loc10 = [None] * (bins-1)
T89_before = [None] * (bins-1)


for i in range(len(epoch)):
    try:
        loc9[i] = t_cis[i]>t_cis[i][index_T89[i]]
        T89_after[i] = h_dens[i][loc9[i]].mean() + o_dens[i][loc9[i]].mean()
    except IndexError:
        print(np.array(loc7[i]).shape)
        print(loc9[i])
        T89_after[i] = [None]
        print(T89_after[i])
for i in range(len(epoch)):
    try:
        loc10[i] = t_cis[i]<t_cis[i][index_T89[i]]
        T89_before[i] = h_dens[i][loc10[i]].mean() + o_dens[i][loc10[i]].mean()
    except IndexError:
        print(np.array(loc10[i]).shape)
        print(loc10[i])
        T89_before[i] = [None]
        print(T89_before[i])

T89_diff = [None] * (bins-1)

for i in range(len(epoch)):
    try:
        T89_diff[i] = abs(T89_after[i]-T89_before[i])
    except TypeError:
        print(T89_after[i])
        print(T89_before[i])
        T89_diff[i] = [None]

for i in range(len(epoch)):
    timediff_T89[i] = timediff_T89[i][0]

x = [0] * (24)
y = [0] * (24)
for i in range(24):
    x[i] = timediff_T89[i].seconds
    y[i] = T89_diff[i]
z = np.poly1d(np.polyfit(x,y,1))
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
r_sq = r_value**2
sig = std_err*math.sqrt(len(x))
sig2 = sig**2
mx = sum(x)/len(x)
sx2 = [None] * (24)
for i in range(24):
    sx2[i] = ((x[i]-mx)**2)
sum_sx2 = sum(sx2)
std_intercept = std_err * math.sqrt(1./len(x) + mx**2/sum_sx2)
sig_intercept = std_intercept * math.sqrt(len(x))
std_slope = std_err * math.sqrt(1./sum_sx2)
sig_slope = std_slope * math.sqrt(len(x))
#z = np.poly1d(np.polyfit(x,y,1))

fig,ax = plt.subplots()
plt.plot(x,y,'bo',x,z(x),'b-')
plt.errorbar(x, y, yerr = std_err, fmt = 'bo')
plt.xlabel('Time Difference(s)', fontsize = 14)
plt.ylabel('Mass Asymmetry ($cm^-$$^3$)')
plt.suptitle('Mass Asymmetry by Time Difference in Plasma Sheet Crossing', fontsize = 16)
plt.show()


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
