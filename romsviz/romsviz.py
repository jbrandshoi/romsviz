
# ============================================================================================
# TODO (enhance): Consider implementing a namelist approach for info such as dim names, titles/labels and colormaps
# TODO (enhance): Add support for creating subplots in existing figure passed to the various plotting methods
# TODO (enhance): Consider supporting dpeth-horizontal slices not parallell to coordinate axes
# TODO (enhance): Support user inputting **plot_kwargs and pass onto the plot/contourf/... functions
# TODO (enhance): Better error handling ensuring correct dimension limits from user before calling get_var()
# TODO (enhance): Support animation of vertical crossections too
# TODO (issue): Specifying range for both x and y is currently "allowed" in depth_csection()
# ============================================================================================

# general modules
import json
import pandas as pd
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1
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
                      "z": ["s_rho", "s_w"],
                      "t": ["ocean_time"]}
        self.cmaps = {"temp": cmocean.cm.thermal, "salt": cmocean.cm.haline,
                      "default": plt.cm.viridis}
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
        ranged_coors = ["t"]
        self.check_range_dims(var, ranged_coors)
        
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
        
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
        fig.tight_layout()
        
        return fig, ax
    
    def depth_time_contour(self, var_name, figax=None, **limits):
        """Function docstring..."""
        var = self.get_var(var_name, **limits)
        ranged_coors = ["z", "t"]
        self.check_range_dims(var, ranged_coors)
        
        z = self.get_sdepths(var)
        
        # plot depth time contour with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        
        try:
            cs = ax.contourf(var.time, z, var.data.transpose(), cmap=self.cmaps[var.name])
        except KeyError:
            cs = ax.contourf(var.time, z, var.data.transpose(), cmap=self.cmaps["default"])
        
        # colorbar stuff
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2.5%", pad=0.15)
        cbar_label = var.attr_to_string(var.meta, "units")
        cb = plt.colorbar(cs, cax=cax, label=cbar_label, orientation="vertical")
        
        # title and labels
        name = var.attr_to_string(var.meta, ["long_name", "standard_name"])
        limits_str = var.lims_to_str(exclude=[var.time_name, "s_rho"])
        ax.set_title("{} {}".format(name, limits_str))
        ax.set_ylabel("Depth [m]")
        ax = self._set_default_txtprop(ax)
        
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")  # rotate xlabel
        fig.tight_layout()
        
        return fig, ax
    
    def vertical_csection(self, var_name, figax=None, **limits):
        """Function docstring..."""
        var = self.get_var(var_name, **limits)
        ranged_coors = ["x", "y", "z"]
        self.check_range_dims(var, ranged_coors)
        
        z = self.get_sdepths(var)  # depth of s levels associated with <var>
        
        # plot time series with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        
        try:
            cs = ax.contourf(var.time, z_rho, var.data.transpose(), cmap=self.cmaps[var.name])
        except KeyError:
            cs = ax.contourf(var.time, z_rho, var.data.transpose(), cmap=self.cmaps["default"])
        
        # colorbar stuff
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(ax)
        cax = divider.append_axes("right", size="2.5%", pad=0.15)
        cbar_label = var.attr_to_string(var.meta, "units")
        cb = plt.colorbar(cs, cax=cax, label=cbar_label, orientation="vertical")
        
        # title and labels
        name = var.attr_to_string(var.meta, ["long_name", "standard_name"])
        limits_str = var.lims_to_str(exclude=[var.time_name, "s_rho"])
        ax.set_title("{} {}".format(name, limits_str))
        ax.set_ylabel("Depth [m]")
        ax = self._set_default_txtprop(ax)
        
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")  # rotate xlabel
        fig.tight_layout()
        
        return fig, ax
    
    def field_2d(self):
        """Function docstring..."""
        raise NotImplementedError("Needs to be overwritten by subclass!")
    
    def images_to_mp4(self, method="convert"):
        """Function docstring..."""
        raise NotImplementedError
    
    def check_range_dims(self, var, answers):
        """Function docstring..."""
        range_dims = self._get_range_dims(var.dim_names, var.lims)
        all_possible = [d for c in answers for d in self.coors[c] if d in var.dim_names]
        
        for r_dim in range_dims:
            if r_dim not in all_possible:
                raise ValueError("{} cannot span multiple indices here!".format(r_dim))
                
    def get_sdepths(self, var):
        """Function docstring..."""
        # get necessary variables for the vertical grid
        h = self.get_var("h").data
        H_c = float(self.get_var("hc").data)
        vtrans = float(self.get_var("Vtransform").data)
        
        # prepare indices for what slice to extract
        x_name = var.identify_dim(self.coors["x"])
        y_name = var.identify_dim(self.coors["y"])
        z_name = var.identify_dim(self.coors["z"])
        xlim = var.extract_lim(x_name)
        ylim = var.extract_lim(y_name)
        zlim = var.extract_lim(z_name)
        slices = self._lims_to_slices([zlim, ylim, xlim])
        
        # variable defined at rho s-levels
        if z_name == "s_rho":
            C = self.get_var("Cs_r").data
            z = roppy.sdepth(h, H_c, C, Vtransform=vtrans)[slices]
        
        # variable defined at w s-levels
        elif z_name == "s_w":
            C_w = self.get_var("Cs_w").data
            z = roppy.sdepth(h, H_c, C_w, Vtransform=vtrans)[slices]  # w s-levels
        
        return z.squeeze()
    
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
        
