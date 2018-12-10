
import sys
import glob
import datetime as dt
import numpy as np
import netCDF4

# Things left to fix (both real issues and enhancements. The number of edge cases is insane!
# You should see the length of this list if I didn't remove stuff along the way.
# ============================================================================================
# TODO (enhance): Allow user to input single point coordinates (instead of tuples) if they want
#       to extract data from a single point in any of the dimsnions
# TODO (issue): Allow user to not specify time dimension and get everything (in time)
# TODO (issue): Some more docstrings
# TODO (enhance): Clean up (enhance) self._compute_time_dist()
# TODO (issue): Pretty sure the if j == 0 in _compute_time_dist() will mess things up when input files only contain 1 time entry
# TODO (premium enhance): Replace self._get_var_{1,2,3,4}d() with fancy indexing, i.e. crete an array (of same shape as the output array) indices gotten from the index limits from user (in progress)
# TODO: (enhance) Consider storing full time array (across files) in __init__ (if exists) and use it for dim lims later
# TODO: (enhance) Fix multi-calculation of var_dim_names
# ============================================================================================

class NetcdfOut(object):
    """Class docstring...
    
    Important Assumptions:
        * If called with wildcard or list of multiple filenames, all data in the different
          files are assumed to be of the same structure (all variables are present in all
          files), only difference being the elements in the time array. Also, automatic
          expansion of the wildcard is assumed to list the files in the correct order.
        * There are only one unlimited dimension and that is time.
    """
    def __init__(self, filepath):
        """
        Constructor function that sets attributes and opens all input files.
        Also extracts the dimensions of the data set for later use.
        
        Args:
            filepath (str/list) : Path/wildcard/list to netcdf data file(s)
        """
        self.filepath = filepath
        self.netcdf_list = self._open_data()
        self.dims = {str(k): v for k, v in self.netcdf_list[0].dimensions.iteritems()}
        self.time_name, self.time_dim = self._get_unlimited_dim()
    
    def _open_data(self):
        """
        Function that parses instance attribute self.filepath and interpretates it
        as either a string (can be wildcard) or list and opens all mathcing netcdf files.
        
        Returns:
            netcdf_list (list(str)) : List of open netcdf file objects
        """
        # depending on type supplied by user, create list fo filepaths
        if type(self.filepath) is str:
            filepath_list = glob.glob(self.filepath)  # expand potential wildcard
        
        elif type(self.filepath) is list:
            filepath_list = self.filepath  # keep user inputed list
        
        else:
            raise TypeError("{} needs to be str or list".format(type(self.filepath)))
        
        if len(filepath_list) == 0:
            raise IOError("Invalid path(s) {}!".format(self.filepath))

        filepath_list = sorted(filepath_list)  # sort in case wildcard or user forgot
        netcdf_list = list()                   # to hold netcdf4 objects (open files)
        
        # loop over each filepath and open netcdf file
        for path in filepath_list:
            try:
                print("Opening file {}".format(path))
                netcdf_list.append(netCDF4.Dataset(path))
            
            except Exception as e:
                sys.exit("Error: Failed opening file {}, {}".format(path, e))
        
        return netcdf_list
    
    def _get_unlimited_dim(self):
        """Function that finds the unlimited dimension in the dataset. Returns
        (None, None) if nonunlimited dimension exists).
        
        Returns:
            dim_name (str)          : Name of unlimited dimension
            dim (netCDF4.Dimension) : Dimension object for unlimtied dim
        """
        for dim_name, dim in self.dims.iteritems():
            if dim.isunlimited():
                return dim_name, dim
        
        return None, None
    
    def get_data(self, var_name, **kwargs):
        """
        Function that supervises the fetching of data from a certain netcdf output
        variable. User may define index limits for all dimensions (or only some of
        them) of the variable if e.g. only parts of the simulation time or domain is
        of interest.
        
        Args:
            var_name (string) : Name of variable to be extracted
        
        Kwargs:
            kwargs (str: tuple) : Lower- and upper index limits for a dimension to
                                  the variable. Min index is 0 and max index is the
                                  size of the particular dimension. Limits for a
                                  dimension defaults to (0, dim.size).
        
        Returns:
            data (numpy.ndarray or float) : Data extracted from the variable within
                                            the index limits specified. Only float
                                            when a scalar is requested, else an ndarray.
        """
        var_meta = self._get_var_meta(var_name)             # netcdf4 variable
        self._verify_kwargs(var_meta, **kwargs)
        lims = self._get_dim_lims(var_meta, **kwargs)
        bounds = list(var_meta.shape)
        print(lims)
        
        # the time dimension may span over multiple files (or just one)
        if self.time_name is not None:
            t_dim_idx = self._get_time_dim_idx(var_meta)
            t_lim = self._get_time_lims(var_meta, kwargs[self.time_name])
            t_bound = self._get_num_time_entries()
            lims = self._update_list_for_time(lims, t_lim, var_meta)
            bounds = self._update_list_for_time(bounds, t_bound, var_meta)
            self._verify_dim_lims(lims, bounds, var_meta)
            use_files, t_dist = self._compute_time_dist(*lims[t_dim_idx])
            data_list = list()
            
            # loop through the all data sets and extract data if inside time limits
            for i, ds in enumerate(self.netcdf_list):
                if use_files[i]:
                    lims = self._update_list_for_time(lims, t_dist[i], var_meta)  # current file time limits
                    array = self._get_var_nd(var_name, lims, ds)                  # actual data
                    data_list.append(array)
            
            # concatenate data from (possibly) multiple files along time axis
            data = np.concatenate(data_list, axis=t_dim_idx)
        
        else:
            data = self._get_var_nd(var_name, lims, self.netcdf_list[0])  # use e.g. zeroth dataset
        
        return data
    
    def _get_var_meta(self, var_name):
        """Function that returns the meta data for requested variable.
        
        Args:
            var_name (string) : Name of variable
        
        Returns:
            var_meta (netCDF4.Variable) : Meta data for the variable
        """
        if not var_name in self.netcdf_list[0].variables:
            raise ValueError("Variable {} not in output data".format(var_name))
        
        else:
            return self.netcdf_list[0].variables[var_name]
    
    def _get_var_dim_names(self, var_meta):
        return [str(s) for s in var_meta.dimensions]
    
    def _get_dim_lims(self, var_meta, **kwargs):
        """
        Function that extracts user provided keyword arguments for
        dimension index limits and constructs a list of with the limits
        being in the exact same order as the actual dimensions variable.
        
        Args:
            var_meta (netCDF4.Variable) : Meta data for Variable of interest
        
        Kwargs:
            kwargs (dict) : See kwargs in self.get_var()
        """
        var_dim_names = self._get_var_dim_names(var_meta)
        idx_lims = list()
        
        # loop over all actual dimensions of the variable
        for vd_name in var_dim_names:
            
            # fill limits if missing kwarg for any dimension
            if vd_name not in kwargs.keys():
                kwargs[vd_name] = (0, self.dims[vd_name].size)  # default to entire range
                
            idx_lims.append(kwargs[vd_name])  # store user inputed limits
        
        return idx_lims
    
    def _verify_kwargs(self, var_meta, **kwargs):
        """
        Function that raises error if not all dimension names
        specified in kwargs are in valid dimensions for the variable.
        
        Kwargs:
            kwargs (dict) : See kwargs in self.get_var()
        """
        var_dim_names = self._get_var_dim_names(var_meta)
        
        for key in kwargs.keys():
            if key not in var_dim_names:
                raise ValueError("Variable {} has no dimension {}!".format(
                                 var_meta.name, key))
        
    def _verify_dim_lims(self, lims, bounds, var_meta):
        """
        Function that verifies if the user provided index limits on
        the requested variable are within bounds, raises error if not.
        
        Args:
            lims (list)                 : List of tuples with index limits
            bounds (list)               : List of upper index bounds
            var_meta (netCDF4.Variable) : Meta data for variable
        """
        var_dim_names = self._get_var_dim_names(var_meta)
        
        # check that specified dimension limits are valid
        for (l_1, l_2), length, v_name in zip(lims, bounds, var_dim_names):
            valid_lims = l_1 >= 0 and l_1 <= length and l_2 >= 0 and l_2 <= length
            
            if not valid_lims:
                raise ValueError("Index limits {} are outside (0, {}) for {}!".format(
                                 (l_1, l_2), length, v_name))
            
            if l_2 < l_1:
                raise ValueError("Lower index {} larger than upper {} for {}!".format(
                                 l_1, l_2, v_name))
        
    def _get_num_time_entries(self):
        """Function that computes total length of time dimension (across files)."""
        return sum(d.variables[self.time_name].shape[0] for d in self.netcdf_list)
    
    def _get_time_dim_idx(self, var_meta):
        """
        Function that finds the index for the time dimension.
        
        Args:
            var_meta (netCDF4.Variable) : Meta data for variable
        """
        var_dim_names = self._get_var_dim_names(var_meta)
        return var_dim_names.index(self.time_name)
        
    def _get_time_lims(self, var_meta, lims):
        """
        Function that handles user provided time limits and returns
        index limits spanning (possibly) over several files.
        
        Args:
            lims (list)                 : Start- and end limits for time (can
                                          be both datetime or indices (int))
            var_meta (netCDF4.Variable) : Meta data for variable
        """
        total_length = self._get_num_time_entries()
        
        if type(lims[0]) is dt.datetime:
            return self._idx_from_dates(*lims)
        
        elif type(lims[0]) is int:
            return lims
        
        else:
            raise TypeError("Invalid type {} for {}".format(type(lims[0]), self.time_name))
    
    def _update_list_for_time(self, list_update, element, var_meta):
        """Function docstring..."""
        var_dim_names = self._get_var_dim_names(var_meta)
        idx = var_dim_names.index(self.time_name)
        list_update[idx] = element
        return list_update

    def _idx_from_dates(self, date_start=None, date_stop=None):
        """Function docstring..."""
        # combine together time arrays from all files
        t_list = [ds.variables[self.time_name] for ds in self.netcdf_list]
        t_raw = np.concatenate([t[:] for t in t_list], axis=0)
        t_dates = netCDF4.num2date(t_raw, t_list[0].units)
        
        def idx_from_date(date_array, date):
            idx = np.where(date_array == date)  # assume only 1 occurence of date
            
            if len(idx[0]) == 0:
                raise ValueError("Date {} not in {}!".format(date, self.time_name))
            
            return idx[0][0]
            
        idx_start = idx_from_date(t_dates, date_start)
        idx_stop = idx_from_date(t_dates, date_stop)
        #t_slice = t_dates[idx_start:idx_stop]
        
        return idx_start, idx_stop
    
    def _compute_time_dist(self, idx_start, idx_stop):
        """Function docstring..."""
        t_per_file = [d.variables[self.time_name].shape[0] for d in self.netcdf_list]
        use_files = [False for _ in t_per_file]  # bool values for relevant files
        idx_total = 0  # to count total indices over all files
        stop_found = False
        
        # compute in what file contains idx_start and idx_stop
        for i in range(len(t_per_file)):
            for j in range(t_per_file[i]):
                if idx_start == idx_total:
                    file_start = i          # file index for time start
                    i_start = j             # start index within that file
                
                if idx_stop == idx_total:
                    if j == 0:
                        file_stop = i - 1         # shift to previous file
                        i_stop = t_per_file[i-1]  # upper limit in prev file
                    
                    else:
                        file_stop = i           # file index for time stop
                        i_stop = j              # stop index within that file
                    
                    stop_found = True
                    break
                
                idx_total += 1
            
            if stop_found:
                break
        
        # special case if final chosen time is overall final time
        if idx_total == sum(t_per_file):
            file_stop = -1
            i_stop = t_per_file[-1]
            
        use_files[file_start] = True
        use_files[file_stop] = True
        # structure to tell, in what file and what index, the specified time starts and stops
        t_dist = [[None, None] for _ in range(len(self.netcdf_list))]  # default (None, None)
        t_dist[file_start][0] = i_start  # file and idx for start time
        t_dist[file_stop][1] = i_stop         # file and idx for stop time
            
        t_dist = tuple(tuple(d) for d in t_dist)  # convert to tuples for safety
        
        # fill use_files with True between "start" and "end" files
        trues = [i for i, el in enumerate(use_files) if el is True]
        
        if len(trues) == 2 and trues[1] - trues[0] != 1:  # done already if not this
            use_files[trues[0]+1:trues[1]] = [True for _ in range(trues[1]-trues[0]-1)]
        
        
        return use_files, t_dist
    
    def _get_var_nd(self, var_name, lims, data_set):
        """Function docstring..."""
        num_dims = len(lims)

        if num_dims == 0:
            # self._get_var_scalar(var_name, data_set)
            array = data_set.variables[var_name].getValue()  # just a scalar
            
        if num_dims == 1:
            array = self._get_var_1d(var_name, lims, data_set)
        
        elif num_dims == 2:
            array = self._get_var_2d(var_name, lims, data_set)
            
        elif num_dims == 3:
            array = self._get_var_3d(var_name, lims, data_set)
        
        elif num_dims == 4:
            array = self._get_var_4d(var_name, lims, data_set)
        
        else:
            raise NotImplementedError("No num_dims > 4 support yet")
        
        return array
        
    def _get_var_1d(self, var_name, lims, data_set):
        """Function docstring..."""
        array = data_set.variables[var_name][lims[0][0]:lims[0][1]]
        return array
        
    def _get_var_2d(self, var_meta, lims, data_set):
        """Function docstring..."""
        array = data_set.variables[var_name][lims[0][0]:lims[0][1],
            lims[1][0]:lims[1][1]]
        return array
    
    def _get_var_3d(self, var_name, lims, data_set):
        """Function docstring..."""
        array = data_set.variables[var_name][lims[0][0]:lims[0][1],
            lims[1][0]:lims[1][1], lims[2][0]:lims[2][1]]
        return array
        
    def _get_var_4d(self, var_name, lims, data_set):
        """Function docstring..."""
        array = data_set.variables[var_name][lims[0][0]:lims[0][1],
            lims[1][0]:lims[1][1], lims[2][0]:lims[2][1], lims[3][0]:lims[3][1]]
        return array
    
    def set_time_name(self, name):
        """Function docstring..."""
        self.time_name = name
    
    def __del__(self):
        """Destructor function closing all files."""
        for data in self.netcdf_list:
            data.close()
        
