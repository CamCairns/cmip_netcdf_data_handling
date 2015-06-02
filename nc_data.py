from netCDF4 import Dataset
from scipy.interpolate import griddata
import numpy as np
import glob
import os
import errno
from netcdftime import utime, datetime

def generate_ensemble_list(experi_list, freq_list, realm_list, vari_list, mount_dir='mountpoint'):
    """ Generates an ensemble list at the model directory level
    
    Args: 
        experi: {AMIP, SPOOKIE}/{experiment}
        freq: frequency (string)
        realm: realm (string)
        vari: variable (string)
        mount: mountpoint directory name
        
    Returns:
        A list of file pathnames
    """
    location = os.path.join('/Users/camcairns/', mount_dir)
    ensemble_list = []
    for experi in experi_list:
        for freq in freq_list:
            for realm in realm_list:
                for vari in vari_list:
                    ensemble_list.append(os.path.join(location, experi, freq, realm, vari))
#                     print ensemble_list
    return ensemble_list

def get_filepath(experi,freq,realm,vari,model,ensemble='r1i1p1',mount_dir='mountpoint',verbose=False):
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
    location = '/Users/camcairns/' + mount_dir + '/'
    files = []
    if ensemble:
        if verbose:
            print "Looking for netcdfs here : \n", os.path.join(location, experi, freq, realm, vari, model, ensemble, "*.nc")
        files = glob.glob(os.path.join(location, experi, freq, realm, vari, model, ensemble, "*.nc"))
    else:
        if verbose:
            print "Looking for netcdfs here: \n", os.path.join(location, experi, freq, realm, vari, model, "*.nc")
        files = glob.glob(os.path.join(location, experi, freq, realm, vari, model, "*.nc"))
    return files

def load_coord_data(files):
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
    nc = Dataset(files[0])
    if nc.variables.has_key('lat'):
        lat= np.squeeze(nc.variables['lat'][:])
    elif nc.variables.has_key('latitude'):
        lat= np.squeeze(nc.variables['latitude'][:])
    else:
        print "No latitude coordinate vector found"
        lat = None

    if nc.variables.has_key('plev'):
        plev= np.squeeze(nc.variables['plev'][:]) # name difference plev, lev
        plev_flag = 1
    elif nc.variables.has_key('lev'):
        plev= np.squeeze(nc.variables['lev'][:])
        plev_flag = 1
    elif nc.variables.has_key('pfull'):
        plev= np.squeeze(nc.variables['pfull'][:])
        plev_flag = 1
    else:
        print "No pressure coordinate vector found"
        plev = None
        plev_flag = 0
    nc.close
#     time= np.squeeze(nc.variables['time'][:])

    return lat, plev, plev_flag

def extract_nc_time(files, model_size):
    """  Gets time, time from a bunch of netcdf files and concatenates them together

    Args:
        files: list of pathways to netcdf variables
        model_size: The timelength of the data, used to preallocate an empty array

    Returns:
        A 1D np array of all the time values concatenated together
    """
    model_size_tkr = 0
    time_vector = np.empty([model_size])*np.nan
    for k in range(len(files)):
        nc = Dataset(files[k],'r')
        if k==0:
            time_units = nc.variables['time'].units
            time_cal = nc.variables['time'].calendar
        time = np.squeeze(nc.variables['time'][:])
        time = time.copy() # make a copy because time is a read only version 
        nc.close
        time_vector[model_size_tkr:model_size_tkr + np.size(time)] = time
        model_size_tkr = model_size_tkr + np.size(time)
    return time_vector, time_units, time_cal

def empty_array_generator(files, time_length, model_dim=None,vari_dim=None,modulo=None):
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
    if modulo:
        padding = (modulo - time_length%modulo)%modulo
    else:
        padding = 0

    if files:
        lat, plev, plev_flag = load_coord_data(files)
        if plev_flag:
            if model_dim:
                nan_array = np.empty([time_length + padding, np.size(plev), np.size(lat), model_dim, vari_dim])*np.nan
            else:
                nan_array = np.empty([time_length + padding, np.size(plev), np.size(lat)])*np.nan
        else: 
            if model_dim:
                nan_array = np.empty([time_length + padding, np.size(lat), model_dim, vari_dim])*np.nan
            else:
                nan_array = np.empty([time_length + padding, np.size(lat)])*np.nan
    else:
        nan_array = None
        print "No file found at that location"
    return nan_array, lat, plev, plev_flag

def extract_nc_data(files, nc_vari, tmp_array, error_limit=1.0e8,zonal_mean=True, verbose=False):
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
    model_size_tkr = 0
    if verbose:
        print 'The model consists of {} .nc file(s)'.format(np.size(files))
    for k in range(len(files)):
        nc = Dataset(files[k],'r')
        temp = np.squeeze(nc.variables[nc_vari][:])
        temp = temp.copy() # make a copy because temp is a read only version 
        nc.close
        if zonal_mean:
            temp[temp>error_limit]= np.nan; # I'd like to make this code genreal for any variable, use dictionary here to assign different tolerances, maybe extract it from the netcdf file?
            temp = np.squeeze(np.nanmean(temp,len(np.shape(temp))-1)); # Take a zonal mean, id like to make this more general as well, automatically pick out the lon dimension?
        tmp_array[model_size_tkr:model_size_tkr+np.size(temp,0),...] = temp
        model_size_tkr = model_size_tkr + np.size(temp,0)
    return tmp_array

def find_model_size(files,nc_variable_name):
    """Finds the total time length of a group of netcdf files

    Args:
        files: list of pathways to netcdf variables
        nc_variable_name: A IPCC CMIP variable shortname (ex. 'ua')
        
    Returns:
        A scalar value, gives the sum of size of the time dimension for each netcdf files in files
    """
    model_size=0;
    for k in range(len(files)):
        nc = Dataset(files[k])
        tmp = np.squeeze(nc.variables[nc_variable_name][:]) # loading dimensions backwards [time, pfull, lat, lon]
        model_size = model_size + np.size(tmp,0);
        nc.close
    return model_size

def interp_data(lat_old, plev_old, plev_flag, lat_new, plev_new, array):
    """interpolates data of the form [time, lat, plev] the function interpolates data array onto a new grid [time, lat_new, plev_new]

The function uses the griddata function, this allows NaNs, (particularly boundary NaNs) to be avoided.

    Args:
        lat_old: vector of the latitude values for the orginal data array
        plev_old: vector of the pressure values for the orginal data array
        lat_new: vector of the latitude values for the new data array
        plev_new: vector of the pressure values for the new data array
        array: original data_array with dimensions [time, lat_old, plev,old]

    Returns:
        A numpy array of interpolated data with dimensions
    """
    if plev_flag:
        array_interp = np.empty([np.size(array,0), len(plev_new), len(lat_new)])
        xold, yold = np.meshgrid(lat_old,plev_old);
        xnew, ynew = np.meshgrid(lat_new,plev_new);
    else:
        array_interp = np.empty([np.size(array,0), len(lat_new)])

    for m in xrange(np.size(array,0)):
        if np.sum(~np.isnan(array[m,...]))==0:  # Some data in the cmip archive is a month of all NaNs, here we set the output as all NaN for that month
            array_interp[m,...] = np.nan
        else:
            mask = ~np.isnan(array[m,...])
            if plev_flag:
                points =  zip(xold[mask], yold[mask])
                temp = array[m,:,:]
                array_interp[m,:,:] = griddata(points, temp[mask],(xnew,ynew),method='cubic')
            else:
                points = lat_old[mask]
                temp = array[m,:]
                array_interp[m,:] = griddata(points, temp[mask],lat_new,method='cubic')

    return array_interp
        
def write_nc(input_lat, input_latb, input_plev, plev_flag, input_array, time_vector, time_units, time_cal, save_path, model_size, experi, freq, realm, vari, model):
    """ Writes the zonal mean data that has been extracted out as a netcdf file, saves in a new directory structure

        files: list of pathways to netcdf variables
        tmp_array: An preallocated numpy vector of appropriate length

    Returns:
        Nothing, a netcdf file is written at the path specified at save_path
    """
    units_dict = {'ua': 'm/s', 'va': 'm/s', 'uas': 'm/s', 'vas': 'm/s', 'ta': 'K', 'tas': 'K', 'hur': '%', 'hus': '1'}
    # Initiate NETCDF4
    rootgrp = Dataset(save_path, 'w', format='NETCDF4')

    # Set Dimensions
    if plev_flag:
        plev = rootgrp.createDimension('plev', len(input_plev))
    else:
        pass
    lat = rootgrp.createDimension('lat', len(input_lat))
    latb = rootgrp.createDimension('latb', len(input_latb))
    time = rootgrp.createDimension('time', None)

    # Set Variables
    times = rootgrp.createVariable('time','f8',('time',))
    if plev_flag:
        plev = rootgrp.createVariable('plev','f8',('plev',))
        tmp = rootgrp.createVariable(vari,'f8',('time','plev','lat'))
    else:
        tmp = rootgrp.createVariable(vari,'f8',('time','lat'))
    latitudes = rootgrp.createVariable('lat','f8',('lat',))
    latitude_bnds = rootgrp.createVariable('lat_bnds','f8',('latb',))
    # two dimensions unlimited.

    # Attributes for nc file
    import time
    rootgrp.description = 'SPOOKIE data with a zonal mean taken and then interpolated on to a uniform lat-lev grid. This nc_files is for MODEL: ' + model + '  FREQ: ' + freq + '  EXPERIMENT: ' + experi + '  REALM: ' + realm + '  VARI: ' + vari
    rootgrp.history = 'Created ' + time.ctime(time.time()) + '. Original data sourced from DKRZ EUCLIPSE sever.'
    rootgrp.source = 'Cameron Cairns; email: cam.cairns1@gmail.com'
    latitudes.units = 'degrees north'
    if plev_flag:
        plev.units = 'hPa'
    else:
        pass
    tmp.units = units_dict[vari]
    if 'since' in time_units:
        times.units = time_units
    else: 
        print "We have an error in times.units" 
    times.calendar = time_cal
    # Write Data
    latitudes[:] = input_lat
    latitude_bnds[:] = input_latb
    times[:] = time_vector
    if plev_flag:
        plev[:] = input_plev
        tmp[0:model_size,0:len(input_plev),0:len(input_lat)] = input_array
    else: 
        tmp[0:model_size,0:len(input_lat)] = input_array
    rootgrp.close()

def mkdir_p(path):
    """Checks to see if a directory path exists
    
    Args:
        path: directory path whose existence needs to be checked.
    """
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def mod_subtract(a,b,N=12):
    """Does subtraction in modulo N (default N=12)
    Used in round_time
    """
    x = (a - b) % N
    return x
    
def round_time(files, model_size,start_month=1):
    """Finds how many months R need to be dropped to initialise each model at start_month
    Function takes in netcdf time data and determines when the model run starts (most start at January (i.e month=1) but some modelling groups start on different months 
    (the UKMO models tend to begin in Dec). It then dtermines the number on months needed to be dropped from the start of the model run to begin the time series at start_month
    
    Args:
        files: directory path to the relevant netcdf data
        model_size: time_length of the model
        start_month (optional): the start month that each model needs to initialise its time series from (1=Jan, 12=Dec), default is January

    Returns:
        R: the remainder, the number of months needed to be dropped (necessarily 0<R<12)
            example: if start_month = 1 and the model originally starts in December then R=1 ==> The first month (Dec) needs to be dropped to start at Jan.
    """
    time_vector, time_units, time_cal = extract_nc_time(files, model_size)
    cdftime = utime(time_units,calendar=time_cal)
    T0 = cdftime.num2date(time_vector[0])
    R = mod_subtract(start_month,T0.month)
    return R

def reshape_data(output_array,plev_flag):
    """Reshapes a data array of form [time, plev OR lat, model] into [months, year, plev OR lat, model] form
    The function assumes the size of output array has been constructed to be divisible by 12 (this is done in round_time and empty_array_generator)

    Args:
        output_array: array to be reshaped
        plev_flag: =1 if field is of form [time plev lat models]; =0 if [time lat models]

    Returns:
        output_array: the reshaped array with shape [months(=12), year, plev OR lat, model]
    """
#     if np.shape(output_array,0)%12 != 0 
        # Write exception here, really should never happen but a useful check
    if plev_flag:
        output_array = np.reshape(output_array, (12,np.size(output_array,0)/12,np.size(output_array,1),np.size(output_array,2),np.size(output_array,3),np.size(output_array,4)),order='F')
    else:
        output_array = np.reshape(output_array, (12,np.size(output_array,0)/12,np.size(output_array,1),np.size(output_array,2),np.size(output_array,3)),order='F')
    return output_array
