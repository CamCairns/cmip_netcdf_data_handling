import numpy as np
import glob
import os
from netCDF4 import Dataset

def get_filepath(experi, mount_dir='mountpoint', verbose=False):
    """Gets a list of the SPOOKIE filepaths from the directory structure I have created,
    Directory structure has form SPOOKIE/experi/freq/realm/vari/model/. Operates using a directory mounted up using FUSE OSX
    
    Args:
        experi: {AMIP, SPOOKIE}/{experiment}
        freq: frequency (string)
        realm: realm (string)
        vari: variable (string)
        model: model (string)
        ensemble: ensemble number (string,optional)
        mount: mountpoint directory name
        
    Returns:
        A list of file pathnames
    """
    location = '/Users/camcairns/' + mount_dir + '/atm_dycores_jucker_strat/'
    dirs = []
    if verbose:
        print "Looking for netcdfs here : \n", os.path.join(location, experi, "run*")
    dirs = glob.glob(os.path.join(location, experi, "run*"))
    model_size = len(dirs)
    return dirs, model_size

def load_coord_data(dirpath,filename):
    """Extracts coordinate data from a bunch of nc files

    Given a list of file pathnames the function takes the first filename and extracts lat and pressure (if it exists). 
    Currently the function looks for pressure levels to be named either 'plev', 'lev' or 'pfull' and if found sets plev_flag = 1. If it does not find any field with this title it 
    assumes the data is a surface field (eg. tas, uas etc.) and the plev_flag = 1. Plev_flag is used in other functions.

    Args:
        files: list of pathways to netcdf variables
    
    Returns
        lat: A vector of the lat values
        plev: A vector of pressure level values (set to None if no pressure dimension exists)
        plev_flag: =1 if plev exists, =0 otherwise
    """
    filepath = os.path.join(dirpath,filename)
    nc = Dataset(filepath,'r')
    lat= np.squeeze(nc.variables['lat'][:])

    if nc.variables.has_key('plev'):
        plev= np.squeeze(nc.variables['plev'][:]) # name difference plev, lev
        plev_flag = 1
    elif nc.variables.has_key('lev'):
        plev= np.squeeze(nc.variables['lev'][:])
        plev_flag = 1
    elif nc.variables.has_key('pfull'):
        plev= np.squeeze(nc.variables['pfull'][:])
        plev_flag = 1
    elif nc.variables.has_key('level'):
        plev= np.squeeze(nc.variables['level'][:])
        plev_flag = 1
    else:
        plev = []
        plev_flag = 0
    nc.close
    
#     if nc.variables.has_key('lon'):
#         lon= np.squeeze(nc.variables['lon'][:])
#         print lon
#     else:
#         lon = []
    return plev, lat, plev_flag

def empty_array_generator(dimension_list):
    """Extracts coord data and preallocates a NaN array of appropriates size

    Given a list of file pathnames a time_length and (optinal) a 4th model dimension the function generates a NaN array of the right size [time_length, plev, lat, model_dim].
    If no model dimension is given the array takes on the shape [time_length, plev, lat]. 
    
    The modulo provision allows an empty array of size >= time_length to be gernerated such that S, the size of the actual array generated has the property S%modulo = 0 
    (ie. it pads the array out by rounding up to the next largest whole multiple of modulo)

    Args:
        files: list of pathways to netcdf variables
        time_length: an integer value for how long the time dim of the array should be
        model_dim (optional): the size of the 4th dimension, if none is provided, no 4th dimension is produced
        modulo (optional): an integer value such that S%modulo = 0, where S is the size of the array generated. If modulo=None no padding is done.
    Returns
        nan_array = a NaN array of appropriate size, used for data preallocation
        lat: A vector of the lat values
        plev: A vector of pressure level values (set to None if no pressure dimension exists)
        plev_flag: =1 if plev exists, =0 otherwise
    """
    # We assume here that the [lat plev] size of all models for a particular variable are equal/consistent!
    dimension_list = filter(None,dimension_list)
    nan_array = np.empty(dimension_list)*np.nan
    return nan_array

def get_ensemble_time_dim(time_length, experi_list, mount_dir='mountpoint', verbose='False'):
    if time_length in ('max','min'):
        model_size_list = []
        for experi in experi_list:
            dirs, model_size = get_filepath(experi, mount_dir=mount_dir)
            if dirs:
                model_size_list.append(model_size)
            else:
                model_size_list.append(np.nan)
        if time_length=='max':
            time_length = max(model_size_list) 
        else:
            time_length = min(model_size_list)
        if verbose:
            print 'Maximum model time_length is {}'.format(max(model_size_list))
            print 'Minimum model time_length is {}'.format(min(model_size_list))
    else:
        pass
    print "The time length of the output array is ", time_length
    return time_length

def extract_nc_data(time_length, experi_list, vari_list, filename, mount_dir='mountpoint', verbose=False):
    """Extracts data from a bunch of nc files

    Given a list of files pathname to a bunch of netcdf files, function, opens each files, places into a preallocated numpy array "tmp_array" into a larger array
    Currently also takes a zonal mean but I would like to make this optional

    Args:
        files: list of pathways to netcdf variables
        vari: A IPCC CMIP variable shortname (ex. 'ua')
        tmp_array: An preallocated numpy array of appropriate size
        error_limit: floating point value, data entries greater then this size will be marked as NaNs
        verbose: Verbosity flag, if true will print out some (hopefully!) helpful statements
    Returns:
        A numpy array of all the netcdf files concatenated together
    """
    time_length = get_ensemble_time_dim(time_length, experi_list, mount_dir=mount_dir)
    dirs, model_size = get_filepath(experi_list[0],mount_dir=mount_dir)
    plev, lat, plev_flag = load_coord_data(dirs[0],filename)
    output_array = empty_array_generator([time_length, len(plev), len(lat), len(experi_list), len(vari_list)])
    for i1, vari in enumerate(vari_list):
        for i2, experi in enumerate(experi_list):
            dirs, model_size = get_filepath(experi, mount_dir=mount_dir)
            run_size = 1 #TODO: Write function to extract this automatically.
            for i3, dir in enumerate(dirs):
                filepath = os.path.join(dir,filename)
                nc = Dataset(filepath,'r')
                tmp = np.squeeze(nc.variables[vari][:])
                tmp = tmp.copy() # make a copy because tmp is a read only version 
                nc.close
                output_array[i3*run_size:(i3+1)*run_size,:,:,i2,i1] = tmp
                nc.close
    return output_array, plev, lat
