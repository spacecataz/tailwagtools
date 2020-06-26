#!/usr/bin/env python

'''
A quick test to see if minute resolution OMNI data is available.
'''

import datetime as dt
import matplotlib.pyplot as plt
from spacepy import omni, time

from spacepy.plot import style

style()

# Create an array of minute-spaced timedeltas.
t = [dt.datetime(2002, 10, 23, 0, 0)+dt.timedelta(minutes=i)
     for i in range(60*12)]

# Convert to ticktocks:
ticks = time.Ticktock(t, 'ISO')

# Fetch hour-level omni:
om1 = omni.get_omni(ticks)

# Fetch minute-level omni:
om2 = omni.get_omni(ticks, dbase='qd1min')

# Compare hour- to minute-resolution data:
plt.plot(om1['UTC'], om1['BzIMF'], label='Hourly data')
plt.plot(om2['UTC'], om2['BzIMF'], label='Minute data')
plt.ylabel('IMF B$_Z$ ($nT$)')
plt.xlabel('Time')
plt.legend(loc='best')
plt.tight_layout()

plt.show()
