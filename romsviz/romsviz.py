
import ncout
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

class RomsViz(ncout.NetcdfOut):
    def __init__(self, filename):
        super(RomsViz, self).__init__(filename)
        plt.style.use("seaborn-deep")
        self.default_title_fn = "serif"
        self.default_label_fn = "serif"
        self.default_title_fs = 20
        self.default_label_fs = 15
    
    def time_series(self, var_name, figax=None, **limits):
        """Function docstring..."""
        var = self.get_var(var_name, **limits)
        
        if len(var.data.shape) > 1:
            raise ValueError("Only {} can be non-single index!".format(self.time_name))
        
        # generate a string from dimension limits
        limits_str = var.generate_lims_string(exclude=[var.time_name])
        
        # plot time series with dates on the x-axis
        fig, ax = self._get_figax(figsize=(12,5), figax=figax)
        ax.grid(True)
        ax.plot(var.time, var.data, linewidth=1)
        ax.set_title("{}{}".format(var.meta.long_name.title(), limits_str))
        ax.set_ylabel("{} [{}]".format(var.name, var.meta.units))
        ax = self._set_default_text(ax)
        fig.autofmt_xdate()
        fig.tight_layout()
        
        return fig, ax
    
    def plot_depth_time_series(self):
        raise NotImplementedError
    
    def images_to_mp4(self, method="convert"):
        raise NotImplementedError
    
    def plot_field2d(self):
        raise NotImplementedError("Needs to be overwritten by subclass!")
    
    def _get_figax(self, figsize=(12,7), figax=None):
        if figax is None:
            return plt.subplots(figsize=figsize, facecolor="white")
        
        else:
            return figax[0], figax[1]
    
    def _set_default_text(self, ax):
        ax.title.set_fontname(self.default_title_fn)
        ax.xaxis.get_label().set_fontname(self.default_label_fn)
        ax.yaxis.get_label().set_fontname(self.default_label_fn)
        ax.title.set_fontsize(self.default_label_fs)
        ax.xaxis.get_label().set_fontsize(self.default_label_fs)
        ax.yaxis.get_label().set_fontsize(self.default_label_fs)
        return ax
    
    def __str__(self):
        raise NotImplementedError("Hax!")
