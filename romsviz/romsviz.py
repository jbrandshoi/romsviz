
import ncout
import pandas as pd
import matplotlib.pyplot as plt
import datetime as dt

class RomsViz(ncout.NetcdfOut):
    def __init__(self, filename):
        super(RomsViz, self).__init__(filename)
        self.default_title_fn = "serif"
        self.default_label_fn = "serif"
        self.default_title_fs = 20
        self.default_label_fs = 15
    
    def plot_time_series(self, var_name, **indices):
        data = self.get_data(var_name, **kwargs).squeeze()
        
        if len(data.shape) > 1:
            raise ValueError("Only {} can be non-single index!".format(self.time_name))
        
        # fetch time array for relevant period
        t_start = indices[self.time_name][0]
        t_stop = indices[self.time_name][1]
        
        if type(t_start) is dt.datetime:
            t_start = self._idx_from_date(t_start)
            t_stop = self._idx_from_date(t_stop)
        
        t = self.time[t_start:t_stop+1]  # time in dates for relevant period
        
        indices_str = ""
        
        for k, v in indices.iteritems():
            if k != self.time_name:
                indices_str += "{}: {}, ".format(k, v)
            
        var_meta = self._get_var_meta(var_name)  # to get long_name and units
        fig, ax = self.get_figax()
        ax.plot(t, data, linewidth=1, color="k")
        ax.grid(True)
        title = ax.set_title("{}, {}".format(var_meta.long_name.title(), indices_str))
        ylabel = ax.set_ylabel("{} [{}]".format(var_name, var_meta.units))
        ax = self.set_default_text(ax)
        fig.autofmt_xdate()
        fig.tight_layout()
        
        return fig, ax
    
    def plot_depth_time_series(self):
        raise NotImplementedError
    
    def images_to_mp4(self, method="convert"):
        raise NotImplementedError
    
    def plot_field2d(self):
        raise NotImplementedError("Needs to be overwritten by subclass!")
    
    def get_figax(self):
        return plt.subplots(figsize=(12,7), facecolor="white")
    
    def set_default_text(self, ax):
        ax.title.set_fontname(self.default_title_fn)
        ax.xaxis.get_label().set_fontname(self.default_label_fn)
        ax.yaxis.get_label().set_fontname(self.default_label_fn)
        ax.title.set_fontsize(self.default_label_fs)
        ax.xaxis.get_label().set_fontsize(self.default_label_fs)
        ax.yaxis.get_label().set_fontsize(self.default_label_fs)
        return ax
    
    def __str__(self):
        raise NotImplementedError("Hax!")
