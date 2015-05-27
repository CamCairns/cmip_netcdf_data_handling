from netCDF4 import Dataset
from scipy.interpolate import griddata
import numpy as np
import glob
import os
import errno

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
    if nc.variables.has_key('latitude'):
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

def extract_nc_data(files, nc_vari, tmp_array, error_limit=1.0e8,zonal_mean=True):
    """Extracts data from a bunch of nc files

    Given a list of files pathname to a bunch of netcdf files, function, opens each files, places into a preallocated numpy array "tmp_array" into a larger array
    Currently also takes a zonal mean but I would like to make this optional

    Args:
        files: list of pathways to netcdf variables
        vari: A IPCC CMIP variable shortname (ex. 'ua')
        tmp_array: An preallocated numpy array of appropriate size
        error_limit: floating point value, data entries greater then this size will be marked as NaNs

    Returns:
        A numpy array of all the netcdf files concatenated together
    """
    model_size_tkr = 0
    print 'The model consists of %d .nc file(s)' % np.size(files)
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
    
def extract_nc_time(files, tmp_array):

    """  Gets time, time from a bunch of netcdf files and concatenates them together

        files: list of pathways to netcdf variables
        tmp_array: An preallocated numpy vector of appropriate length

    Returns:
        A numpy array of all the netcdf files concatenated together
    """
    model_size = 0
    for k in range(len(files)):
        print 'Extracting file number %d of %d' % (k+1, np.size(files))
        nc = Dataset(files[k])
        temp = np.squeeze(nc.variables['time'][:])
        temp_copy = temp.copy()
        nc.close
        tmp_array[model_size:model_size+np.size(temp)] = temp_copy
        model_size = model_size + np.size(temp_copy)
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
    
def get_filepath(category,experi,freq,realm,vari,model,ensemble='r1i1p1',mount_dir='mountpoint'):
    """Gets a list of the SPOOKIE filepaths from the directory structure I have created,
    
    Directory structure has form SPOOKIE/experi/freq/realm/vari/model/. Operates using a directory mounted up using FUSE OSX
    
    Args:
        category: AMIP, SPOOKIE
        experi: experiment (string)
        freq: frequency (string)
        realm: realm (string)
        vari: variable (string)
        model: model (string)
        ensemble: ensemble number (string,optional)
        mount: mountpoint directory name
        
    Returns:
        A list of file pathnames
    """
    location = '/Users/camcairns/' + mount_dir + '/' + category + '/'
    files = []
    if ensemble:
        print "Looking for netcdfs here: \n", location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/' + str(ensemble) + '/*.nc'
        files = glob.glob(location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/' + str(ensemble) + '/*.nc')
    else:
        print "Looking for netcdfs here: \n", location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/*.nc'
        files = glob.glob(location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/*.nc')

    
    return files

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

    Raises:
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
        
def write_nc(input_lat, input_latb, input_plev, plev_flag, input_array,  time_array, time_units, time_cal, save_path, model_size, experi, freq, realm, vari, model):
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
        plev = rootgrp.createDimension('plev', 17)
    else:
        pass
    lat = rootgrp.createDimension('lat', 64)
    latb = rootgrp.createDimension('latb', 65)
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
    times.units = time_units
    times.calendar = time_cal
    
    # Write Data
    latitudes[:] = input_lat
    latitude_bnds[:] = input_latb
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
