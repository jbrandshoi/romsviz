
import glob
import datetime as dt
import matplotlib.pyplot as plt
import romsviz.ncout as no

filename = "/lustre/storeB/users/josteinb/metroms_run/barents-2.5km/tmp_3M_data/ocean_avg_000*.nc"
d0 = dt.datetime(2017,10,30,12)
d1 = dt.datetime(2017,10,30,12)
limits = dict(xi_rho=(50,65), eta_rho=(540,545), ocean_time=(d0,d1))

ncout = no.NetcdfOut(filename)
#var = ncout.get_data("zeta", **limits)
time = ncout.get_data("ocean_time", ocean_time=(d0,d1))

for t in time:
    date = dt.datetime(1970,1,1) + dt.timedelta(seconds=t)
    print(date)

#print(var.shape)
