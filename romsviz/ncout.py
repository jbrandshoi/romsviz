
# ============================================================================================
# TODO: (issue) Some more docstrings
# TODO: (enhance) Consider storing full time array (across files) in __init__ (if exists) and use it for dim lims later
# TODO: (enhance) Possibly move datatype checking in _get_dim_lims() to _verify_kwargs()
# TODO: (issue) Using netCDF4.num2date() with var.units is risky
# ============================================================================================

import sys
import glob
import logging
import collections
import datetime as dt
import numpy as np
import netCDF4
import outvar

class NetcdfOut(object):
    """Class docstring...
    
    Important Assumptions:
        * If called with wildcard or list of multiple filenames, all data in the different
          files are assumed to be of the same structure (all variables are present with same
          shape in all files), only difference being the elements in the time array. Also,
          automatic expansion of the wildcard is assumed to list the files in the correct order.
        * There are only one unlimited dimension and that is time.
    """
    def __init__(self, filepath, debug=False):
        """
        Constructor function that sets attributes and opens all input files.
        Also extracts the dimensions of the dataset for later use. Also
        enables/disables logging for debugging purposes.
        
        Args:
            filepath (str/list) : Path/wildcard/list to netcdf data file(s)
            debug (bool)        : True/False to turn on/off debug mode
        """
        if debug:
            logging.basicConfig(level=logging.NOTSET)
            logging.info("debug mode enabled!")
        
        else:
            logging.basicConfig(level=logging.CRITICAL)
            
        self.filepath = filepath
        self.filepaths, self.netcdfs = self.open_data()
        self.dims = {k.encode("utf8"): v for k, v in self.netcdfs[0].dimensions.items()}
        self.default_lim = (None, None)
        self.time_name, self.time_dim = self._get_unlimited_dim()
            
    def open_data(self):
        """
        Function that parses instance attribute self.filepath and interpretates it
        as either a string (can be wildcard) or list and opens all mathcing netcdf files.
        
        Returns:
            netcdf_list (list(str)) : List of open netcdf file objects
        """
        # depending on type supplied by user, create list fo filepaths
        if type(self.filepath) is str:
            filepaths = glob.glob(self.filepath)  # expand potential wildcard
        
        elif type(self.filepath) is list:
            filepaths = self.filepath             # keep user inputed list
        
        else:
            raise TypeError("{} must be str or list".format(self.filepath))
        
        if len(filepaths) == 0:
            raise IOError("Invalid path(s) {}!".format(self.filepath))

        filepaths = sorted(filepaths)  # sort in case wildcard or user forgot
        netcdfs = list()               # to hold netcdf4 objects (open files)
        
        # loop over each filepath and open netcdf file
        for path in filepaths:
            try:
                logging.debug("opening file {}".format(path))
                netcdfs.append(netCDF4.Dataset(path, "r"))
            
            except Exception as e:
                sys.exit("Error: Failed opening file {}, {}".format(path, e))
        
        return filepaths, netcdfs
    
    def _get_unlimited_dim(self):
        """Function that finds the unlimited dimension in the dataset. Returns
        (None, None) if no nunlimited dimension is found in the dataset).
        
        Returns:
            dim_name (str)          : Name of unlimited dimension
            dim (netCDF4.Dimension) : Dimension object for unlimtied dim
        """
        for dim_name, dim in self.dims.items():
            if dim.isunlimited():
                return dim_name, dim
        
        return None, None
    
    def get_var(self, var_name, **limits):
        """
        Function that supervises the fetching of data from a certain netcdf output
        variable. User may define index limits for all dimensions (or only some of
        them) of the variable if e.g. only parts of the simulation time or domain is
        of interest.
        
        Args:
            var_name (string) : Name of variable to be extracted
        
        Kwargs:
            limits (str: tuple) : Lower- and upper index limits for a dimension to
                                  the variable. Min index is 0 and max index is the
                                  size of the particular dimension. Limits for a
                                  dimension defaults to (0, dim.size).
        
        Returns:
            var (outvar.OutVar) : Custom variable object containing info about the
                                  extracted variable including a data array within
                                  the index limits specified. If the variable has a
                                  time dimension, the object contains a datetime array
                                  for the relevant time range of the extracted variable.
        """
        logging.debug("extracting variable {}".format(var_name))
        logging.debug("user supllied dimension limits: {}".format(limits))
        
        # store info in OutVar object and verify user inputed dimension limits
        var = outvar.OutVar()
        var.name = var_name
        var.meta = self._get_var_meta(var.name)
        var.dim_names = [s.encode("utf8") for s in var.meta.dimensions]
        self._verify_kwargs(var.name, var.dim_names, **limits)
        var.lims = self._get_dim_lims(var.dim_names, **limits)
        var.bounds = list(var.meta.shape)
        var.time_name = self._unlimdim_in_dims(var.dim_names)
        
        # the time dimension may span over multiple files
        if var.time_name is not None:
            self.set_time_array()
            i_td = self._get_time_dim_idx(var.dim_names)
            var.bounds[i_td] = self._get_num_time_entries()
            var.lims[i_td] = self._get_time_lims(var.lims[i_td], var.bounds[i_td])
            self._verify_lims(var.lims, var.bounds, var.dim_names)
            var.use_files, var.t_dist = self._compute_time_dist(*var.lims[i_td])
            data_list = list()
            lims = var.lims[:]
            
            # loop through the all data sets and extract data if inside time limits
            for i, ds in enumerate(self.netcdfs):
                if var.use_files[i]:
                    logging.debug("getting data from file {}".format(self.filepaths[i]))
                    lims[i_td] = var.t_dist[i]
                    slices = self._lims_to_slices(lims)
                    array = self._get_var_nd(var.name, slices, ds)
                    data_list.append(array)
            
            # concatenate data from (possibly) multiple files along time axis
            data = np.concatenate(data_list, axis=i_td)
            var.time = self.time[self._lims_to_slices([var.lims[i_td]])]  # include time slice
        
        # very simple if there's no time dimension
        else:
            self._verify_lims(var.lims, var.bounds, var.dim_names)
            slices = self._lims_to_slices(var.lims)
            data = self._get_var_nd(var.name, slices, self.netcdfs[0])  # use e.g. zeroth dataset

        var.data = data.squeeze()  # finally store the main array in var object
        return var
        
    def _get_var_meta(self, var_name):
        """
        Function that returns the meta data for requested variable.
        
        Args:
            var_name (string) : Name of variable
        
        Returns:
            var_meta (netCDF4.Variable) : Meta data for the variable
        """
        if not var_name in self.netcdfs[0].variables:
            raise ValueError("Variable {} not in output data".format(var_name))
        
        else:
            return self.netcdfs[0].variables[var_name]
    
    def _unlimdim_in_dims(self, vd_names):
        if self.time_name is not None:
            for vd_name in vd_names:
                if vd_name == self.time_name:
                    return vd_name
        
        return None
        
    def _get_dim_lims(self, vd_names, **limits):
        """
        Function that extracts user provided keyword arguments for
        dimension index limits and constructs a list of with the limits
        being in the exact same order as the actual dimensions variable.
        
        Args:
            vd_names (list) : List of actual variable dimension names
        
        Kwargs:
            limits (dict) : See limits in self.get_var()
        """
        idx_lims = list()
        
        # loop over all actual dimensions of the variable
        for vd_name in vd_names:
            
            # fill limits if missing kwarg for any dimension
            if vd_name not in limits.keys():
                lim = self.default_lim  # default to entire range
            
            # if user gives e.g. int, float, dt.datetime
            elif type(limits[vd_name]) not in [tuple, list]:
                lim = (limits[vd_name], limits[vd_name])
            
            # if user has given tuple or list of limits
            else:
                lim = limits[vd_name]
                
            idx_lims.append(lim)  # store user inputed limits
        
        return idx_lims
    
    def _get_range_dims(self, dim_names, lims):
        """Function docstring..."""
        range_dims = list()
        
        for d_name, lim in zip(dim_names, lims):
            if lim[1] != lim[0] or None in lim:
                range_dims.append(d_name)
        
        return range_dims
    
    def _verify_kwargs(self, var_name, vd_names, **limits):
        """
        Function that raises error if not all dimension names
        specified in limits are in valid dimensions for the variable.
        
        Kwargs:
            limits (dict) : See limits in self.get_var()
        """
        for key in limits.keys():
            if key not in vd_names:
                raise ValueError("Variable {} has no dimension {}!".format(var_name, key))
        
    def _verify_lims(self, lims, bounds, vd_names):
        """
        Function that verifies if the user provided index limits on
        the requested variable are within bounds, raises error if not.
        
        Args:
            lims (list)     : List of tuples with index limits
            bounds (list)   : List of upper index bounds
            vd_names (list) : List of actual variable dimension names
        """
        # check that specified dimension limits are valid
        for (l_1, l_2), length, vd_name in zip(lims, bounds, vd_names):
            if (l_1, l_2) == self.default_lim:
                continue  # nothing to check if default limits
                
            valid_lims = l_1 >= 0 and l_1 <= length and l_2 >= 0 and l_2 <= length
            
            if not valid_lims:
                raise ValueError("Index limits {} are outside (0, {}) for {}!".format(
                                 (l_1, l_2), length, vd_name))
                                 
            if l_2 < l_1:
                raise ValueError("Lower index {} larger than upper {} for {}!".format(
                                 l_1, l_2, vd_name))
        
    def _get_num_time_entries(self):
        """Function that counts number of total time entries (across files).
        
        Returns:
            num_time_entries (int) : Number of time elements across all input files
        """
        return sum(d.variables[self.time_name].shape[0] for d in self.netcdfs)
    
    def set_time_array(self):
        """Function that stitches together the time array over all files.
        
        Returns:
            t_dates (np.ndarray (1D)) : Full time array across all files
        """
        t_list = [ds.variables[self.time_name] for ds in self.netcdfs]
        t_raw = np.concatenate([t[:] for t in t_list], axis=0)
        t_dates = netCDF4.num2date(t_raw, t_list[0].units)
        self.time = t_dates
    
    def _get_time_dim_idx(self, vd_names):
        """
        Function that finds the index for the time dimension.
        
        Args:
            vd_names (list) : List of actual variable dimension names
        """
        return vd_names.index(self.time_name)
        
    def _get_time_lims(self, t_lim, total_length):
        """
        Function that handles user provided time limits and returns
        index limits spanning (possibly) over several files.
        
        Args:
            lims (list) : Start- and end limits for time (can
                          be both datetime or indices (int))
        """
        int_types = [int, np.int8, np.int16, np.int32, np.int64]
        
        if t_lim == self.default_lim:
            return (0, total_length - 1)
        
        elif type(t_lim[0]) in int_types or type(t_lim[1]) in int_types:
            return t_lim
            
        elif type(t_lim[0]) is dt.datetime:
            idx_start = self._idx_from_date(t_lim[0])
            idx_stop = self._idx_from_date(t_lim[1])
            return (idx_start, idx_stop)
        
        else:
            raise TypeError("Invalid limits {} for {} (use int/None/datetime)".format(
                t_lim, self.time_name))
    
    def _idx_from_date(self, date):
        """Function docstring..."""
        idx = np.where(self.time == date)  # assume only 1 occurrence of date
        
        if len(idx[0]) == 0:
            raise ValueError("Date {} not in {}!".format(date, self.time_name))
        
        return idx[0][0]
    
    def _compute_time_dist(self, idx_start, idx_stop):
        """Function docstring..."""
        t_per_file = [d.variables[self.time_name].shape[0] for d in self.netcdfs]
        use_files = [False for _ in t_per_file]  # bool values for relevant files
        idx_total = 0       # to count total indices over all files
        
        # compute in what file contains idx_start and idx_stop
        for i in range(len(t_per_file)):
            for j in range(t_per_file[i]):
                if idx_start == idx_total:
                    file_start = i          # file index for time start
                    i_start = j             # start index within that file
                
                if idx_stop == idx_total:
                    file_stop = i           # file index for time stop
                    i_stop = j              # stop index within that file
                    
                idx_total += 1
        
        use_files[file_start] = True
        use_files[file_stop] = True
        
        # structure to tell, in what file and what index, the specified time starts and stops
        t_dist = [[None, None] for _ in range(len(self.netcdfs))]  # default (None, None)
        t_dist[file_start][0] = i_start  # file and idx for start time
        t_dist[file_stop][1] = i_stop    # file and idx for stop time
            
        t_dist = tuple(tuple(d) for d in t_dist)  # convert to tuples for safety
        
        # fill use_files with True between "start" and "end" files
        trues = [i for i, el in enumerate(use_files) if el is True]
        
        if len(trues) == 2 and trues[1] - trues[0] != 1:  # done already if not this
            use_files[trues[0]+1:trues[1]] = [True for _ in range(trues[1]-trues[0]-1)]
        
        return use_files, t_dist
    
    def _lims_to_slices(self, lims):
        """Function docstring..."""
        slices = list()
        
        for l in lims:
            if l[1] is not None:
                slices.append(slice(l[0], l[1] + 1))  # include end point too
            
            else:
                slices.append(slice(l[0], l[1]))
        
        return tuple(slices)
        
    def _get_var_nd(self, var_name, slices, dataset):
        """Function docstring..."""
        return dataset.variables[var_name][slices]
    
    def __del__(self):
        """Destructor function closing all files."""
        for data in self.netcdfs:
            data.close()
        
