
import glob
import datetime as dt
import matplotlib.pyplot as plt
import romsviz.ncout as no
import romsviz.romsviz as rv
"""
filename = "/lustre/storeB/users/josteinb/metroms_run/barents-2.5km/tmp_3M_data/ocean_avg_000*.nc"
#filename = "/lustre/storeB/project/fou/hi/new_norkyst/his/ocean_his.an.201802*.nc"
d0 = dt.datetime(2017,10,1,12)
d1 = dt.datetime(2017,11,10,12)
limits = dict(xi_rho=(50,65), eta_rho=(10,20), ocean_time=(d0,d1))

ncout = no.NetcdfOut(filename)
var = ncout.get_var("zeta", **limits)
#time = ncout.get_var("ocean_time")

for t in var.time:
    print(t)
    
print(var.data.shape)
"""

# test romsviz
filename = "/lustre/storeB/users/josteinb/metroms_run/barents-2.5km/tmp_3M_data/ocean_avg_000*.nc"
d0 = dt.datetime(2017,10,1,12)
d1 = dt.datetime(2017,12,31,12)

rviz = rv.RomsViz(filename)
fig, ax = rviz.depth_time_series("temp", ocean_time=(d0,d1), xi_rho=200, eta_rho=200)
plt.show()

