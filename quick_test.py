#!/usr/bin/env python
from scipy.io.matlab import loadmat
from libsrvf import ndPoint, ndPointset, functions_to_srvfs, calculate_karcher_mean
import numpy as np
import matplotlib.pyplot as plt

m = loadmat("/home/matt/simu_data.mat")
# Specify the test data
fake_data = np.ascontiguousarray(m['f'].T)
sample_times = m['time']

print "pt = ndPoint()"
pt = ndPoint()
print pt.dim()
print pt.norm()

print "ps = Pointset(np.arange(15)"
ps = ndPointset(np.arange(5,dtype=np.float))
print "ps.get_data()"
print ps.get_data()
print "ps.scale(10.0)"
ps.scale(10.0)
print ps.get_data()

print "ps.dim()", ps.dim()
print "ps.npts()", ps.npts()

srvf_version = functions_to_srvfs(fake_data)

fig, axes = plt.subplots(nrows=2)

for orig_func in fake_data:
    axes[0].plot(orig_func,'.-')
axes[0].set_title("Original Functions")
    
for srvf_func in srvf_version:
    axes[1].plot(srvf_func, '.-')
axes[1].set_title("Srvf Functions")

fig2, axes2 = plt.subplots()
karcher_data, karcher_time = calculate_karcher_mean(fake_data)
axes2.plot(karcher_time, karcher_data,'.-')
plt.show()
