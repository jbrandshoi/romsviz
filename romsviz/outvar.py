
class OutVar(object):
    def __init__(self):
        self.name = None
        self.meta = None
        self.lims = None
        self.time_dist = None
        self.use_files = None
        self.dim_names = None
        self.data = None
    
    def __str__(self):
        if self.var_meta is None:
            return self.var_name
        
        else:
            return self.var_meta
