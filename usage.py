
import glob
import datetime as dt
import matplotlib.pyplot as plt
import romsviz.ncout as no
import romsviz.romsviz as rv


filename = "/lustre/storeB/users/josteinb/metroms_run/barents-2.5km/tmp_3M_data/ocean_avg_000*.nc"
#filename = "/lustre/storeB/project/fou/hi/new_norkyst/his/ocean_his.an.201802*.nc"
d0 = dt.datetime(2017,10,1,12)
d1 = dt.datetime(2017,11,10,12)
limits = dict(xi_rho=(50,65), eta_rho=(10,20), ocean_time=(d0,d1))

ncout = no.NetcdfOut(filename, debug=True)
var = ncout.get_var("zeta", **limits)
#time = ncout.get_var("ocean_time")

"""
# test romsviz # 542, 211
filename = "/lustre/storeB/users/josteinb/metroms_run/barents-2.5km/tmp_3M_data/ocean_avg_000*.nc"
grid_filename = "/home/josteinb/metroms_apps/barents-2.5km/grid/barents_grd.nc"
d0 = dt.datetime(2017,10,1,12)
d1 = dt.datetime(2017,10,10,12)

rviz = rv.RomsViz(filename)
#fig, ax = rviz.depth_time_contour("temp", ocean_time=(d0,d1), xi_rho=200, eta_rho=200)
rviz.set_gridfile(grid_filename)
fig, ax = rviz.csection("temp", ocean_time=d0, eta_rho=200)
plt.show()
"""

