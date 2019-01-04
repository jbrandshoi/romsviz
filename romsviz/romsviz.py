
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
        self.x_id = "x"
        self.y_id = "y"
        self.z_id = "z"
        self.t_id = "t"
        self.coors = {self.x_id: ["xi_rho", "xi_u", "xi_v", "xi_psi"],
                      self.y_id: ["eta_rho", "eta_u", "eta_v", "eta_psi"],
                      self.z_id: ["s_rho", "s_w"],
                      self.t_id: ["ocean_time"]}
        
        self.cmaps = {"temp": cmocean.cm.thermal, "salt": cmocean.cm.haline,
                      "default": plt.cm.viridis}  # perhaps replace with json input in future
        self.default_title_fs = 20
        self.default_label_fs = 15
        plt.style.use("seaborn-deep")
        plt.rc("font", family="serif")
    
    def set_gridfile(self, filename):
        self.gridfile = ncout.NetcdfOut(filename)
    
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
    
    def csection(self, var_name, figax=None, **limits):
        """Function docstring..."""
        var = self.get_var(var_name, **limits)
        ranged_coors = ["x", "y", "z", "t"]
        range_dims = self.check_range_dims(var, ranged_coors, upper_limit=2)
        
        any_in = lambda a, b: any(i in b for i in a)
        
        if any_in(range_dims, self.coors[self.x_id]):
            if not hasattr(self, "gridfile"):
                raise ValueError("No grid file!")
            
            x_name = var.identify_dim(self.coors[self.x_id])
            x_lim = var.extract_lim(x_name)
            kw = {x_name: x_lim}
            x_axis = self.gridfile.get_var(x_name, **kw).data
        
        if any_in(range_dims, self.coors[self.y_id]):
            raise ValueError("No grid file!")
        
        if any_in(range_dims, self.coors[self.t_id]):
            x_axis = var.time
            
        if any_in(range_dims, self.coors[self.z_id]):
            y_axis = self.get_sdepths(var)  # depth of s levels associated with <var>
            import numpy as np
            y_axis = y_axis[:,np.where(y_axis==np.min(y_axis))[1][0]]
        
        print(x_axis.shape, y_axis.shape, var.data.shape)
        
        # plot time series with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        
        try:
            print(var.data.fill_value)
            print(np.ma.masked_where(var.data==var.data.fill_value, var.data))
            cs = ax.contourf(x_axis, y_axis, var.data, cmap=self.cmaps[var.name])
        except KeyError:
            cs = ax.contourf(x_axis, y_axis, var.data.transpose(), cmap=self.cmaps["default"])
        
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
    
    def horizontal_csection(self, var_name, figax=None, **limits):
        """Function docstring..."""
        limits["s_rho"] = 41  # default to surface
        raise NotImplementedError()
    
    def images_to_mp4(self, method="convert"):
        """Function docstring..."""
        if method == "convert":
            subprocess.call("convert -loop 0 -delay 10 -hax 1 {} {}".format(wildcard, output_fn))
        
        elif method == "builtin":
            anim = animation.FuncAnimation()
            anim.save(output_fn)
            
        raise NotImplementedError
    
    def check_range_dims(self, var, answers, upper_limit=1):
        """Function docstring..."""
        range_dims = self._get_range_dims(var.dim_names, var.lims)
        allowed = [d for c in answers for d in self.coors[c] if d in var.dim_names]
        print(range_dims)
        print(allowed)
        
        for r_dim in range_dims:
            if r_dim not in allowed:
                raise ValueError("{} cannot span multiple indices!".format(r_dim))
        
        if len(range_dims) > upper_limit:
            raise ValueError("Only {} dimensions can span multiple indices!".format(upper_limit))
            
        return range_dims
                
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
            z = roppy.sdepth(h, H_c, C, Vtransform=vtrans, stagger="rho")[slices]
        
        # variable defined at w s-levels
        elif z_name == "s_w":
            C_w = self.get_var("Cs_w").data
            z = roppy.sdepth(h, H_c, C_w, Vtransform=vtrans, stagger="w")[slices]
        
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
        
