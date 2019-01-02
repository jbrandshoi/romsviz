
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
    
    def extract_lim(self, dim_name):
        """Function docstring..."""
        for i in range(len(self.dim_names)):
            if self.dim_names[i] == dim_name:
                return self.lims[i]
        
        raise ValueError("{} is not a dimension of {}!".format(
                         dim_name, self.var_name))
    
    def identify_dim(self, suggestions):
        for dim_name in self.dim_names:
            for d in suggestions:
                if dim_name == d:
                    return dim_name
        
        raise ValueError("No dimension of {} are in <suggestions>".format(self.name))
    
    def attr_to_string(self, obj, attr):
        """Function docstring..."""
        if type(attr) is str:
            attr = [attr]  # need to be list below
        
        for a in attr:
            val = getattr(obj, a, None)
            
            if val:
                return val.encode("utf8").capitalize()
        
        return "N/A"
            
    def lims_to_str(self, exclude=list()):
        """Function docstring..."""
        lims_str = "("
        
        for d_name, lim in zip(self.dim_names, self.lims):
            if d_name not in exclude:
                if lim[0] == lim[1]:
                    lims_str += "{}: {}".format(d_name, lim[0])
                
                else:
                    lims_str += "{}: {}".format(d_name, lim)
                
                if d_name != self.dim_names[-1]:
                    lims_str += ", "
        
        return lims_str + ")"
    
    def __str__(self):
        if self.var_meta is None:
            return self.var_name
        
        else:
            return self.var_meta

