
# general modules
import json
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

# other ROMS specific modules
import roppy    # needed to convert from sigma- to z-coordinates
import cmocean  # for colormaps

# module(s) part of this package
import ncout  # extracting data from netcdfs

class RomsViz(ncout.NetcdfOut):
    def __init__(self, filename, infofile="romsviz/namelist.json"):
        super(RomsViz, self).__init__(filename)
        self.coors = {"x": ["xi_rho", "xi_u", "xi_v", "xi_psi"],
                      "y": ["eta_rho", "eta_u", "eta_v", "eta_psi"],
                      "z": ["z_rho", "z_w"],
                      "t": ["ocean_time"]}
        plt.style.use("seaborn-deep")
        plt.rc("font", family="serif")
        self.default_title_fs = 20
        self.default_label_fs = 15
    
    def parse_namelist(self, infofile):
        """Function docstring..."""
        with open(infofile, "r") as nl:
            data = json.load(nl)
        
        return data
    
    def time_series(self, var_name, figax=None, **limits):
        """Function docstring..."""
        var = self.get_var(var_name, **limits)
        
        if len(var.data.shape) > 1:
            raise ValueError("Only {} can be non-single index!".format(self.time_name))
        
        # plot time series with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        ax.plot(var.time, var.data, linewidth=1)
        ax.grid(True)
        
        name = var.attr_to_string(var.meta, ["long_name", "standard_name"])
        ylabel = var.attr_to_string(var.meta, "units")
        limits_str = var.lims_to_str(exclude=[var.time_name])
        ax.set_title("{} {}".format(name, limits_str))
        ax.set_ylabel("{}".format(ylabel))
        ax = self._set_default_txtprop(ax)
        
        fig.autofmt_xdate()
        fig.tight_layout()
        
        return fig, ax
    
    def depth_time_series(self, var_name, figax=None, **limits):
        """Function docstring..."""
        # 542, 211
        var = self.get_var(var_name, **limits)
        
        if len(var.data.shape) > 2:
            raise ValueError("Only {} and {} can be non-single index!".format(
                var.time_name, "s_rho"))
        
        # get depth of s levels using roppy
        h = self.get_var("h").data
        H_c = float(self.get_var("hc").data)
        vtrans = float(self.get_var("Vtransform").data)
        x_name = var.identify_dim(self.coors["x"])
        y_name = var.identify_dim(self.coors["y"])
        z_name = var.identify_dim(self.coors["z"])
        y_idx = var.idx_from_lims(x_name)[0]
        x_idx = var.idx_from_lims(y_name)[0]
        
        if z_name == "z_rho":
            C = self.get_var("Cs_r").data
            z_rho = roppy.sdepth(h, H_c, C, Vtransform=vtrans)[:,y_idx,x_idx]  # rho s-levels
        
        elif z_name == "z_w":
            C_w = self.get_var("Cs_w").data
            z_w = roppy.sdepth(h, H_c, C_w, Vtransform=vtrans)[:,y_idx,x_idx]  # w s-levels
        
        
        # plot time series with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        ax.contourf(var.time, z_rho, var.data.transpose(), cmap=cmocean.cm.thermal)
        
        name = var.attr_to_string(var.meta, ["long_name", "standard_name"])
        limits_str = var.lims_to_str(exclude=[var.time_name, "s_rho"])
        ax.set_title("{} {}".format(name, limits_str))
        ax.set_ylabel("Depth [m]")
        ax = self._set_default_txtprop(ax)
        
        fig.autofmt_xdate()
        fig.tight_layout()
        
        return fig, ax
    
    def images_to_mp4(self, method="convert"):
        raise NotImplementedError
    
    def field_2d(self):
        raise NotImplementedError("Needs to be overwritten by subclass!")
    
    def _get_figax(self, figsize=(12,7), figax=None):
        if figax is None:
            return plt.subplots(figsize=figsize, facecolor="white")
        
        else:
            return figax[0], figax[1]
    
    def _set_default_txtprop(self, ax):
        #ax.title.set_fontname(self.default_title_fn)
        #ax.xaxis.get_label().set_fontname(self.default_label_fn)
        #ax.yaxis.get_label().set_fontname(self.default_label_fn)
        ax.title.set_fontsize(self.default_title_fs)
        ax.xaxis.get_label().set_fontsize(self.default_label_fs)
        ax.yaxis.get_label().set_fontsize(self.default_label_fs)
        return ax
    
    def __str__(self):
        raise NotImplementedError("Hax!")
