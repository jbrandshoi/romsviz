
class OutVar(object):
    """Class docstring..."""
    def __init__(self):
        """Constructor docstring..."""
        self.name = None
        self.meta = None
        self.lims = None
        self.time_dist = None
        self.use_files = None
        self.dim_names = None
        self.data = None
        self.time_name = None
        self.time = None
    
    def generate_lims_string(self, exclude=list()):
        """Function docstring..."""
        lims_str = ""
        
        for d_name, lim in zip(self.dim_names, self.lims):
            if d_name not in exclude:
                if lim[0] == lim[1]:
                    lims_str += ", {}: {}".format(d_name, lim[0])
                
                else:
                    lims_str += ", {}: {}".format(d_name, lim)
        
        return lims_str
    
    def __str__(self):
        if self.var_meta is None:
            return self.var_name
        
        else:
            return self.var_meta
