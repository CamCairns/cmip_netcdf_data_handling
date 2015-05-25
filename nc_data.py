from scipy.io import netcdf_file
from netCDF4 import Dataset
from scipy.interpolate import griddata
import numpy as np
import glob
import os
import errno

def load_coord_data(files):
    nc = netcdf_file(files[0])
    if nc.variables.has_key('lat'):
        lat= np.squeeze(nc.variables['lat'][:])

    if nc.variables.has_key('plev'):
        plev= np.squeeze(nc.variables['plev'][:]) # name difference plev, lev
    elif nc.variables.has_key('lev'):
        plev= np.squeeze(nc.variables['lev'][:])
    elif nc.variables.has_key('pfull'):
        plev= np.squeeze(nc.variables['pfull'][:])
    else:
        print "treating variable " + vari + " as a surface field"
        plev = None
    nc.close
#     time= np.squeeze(nc.variables['time'][:])
    return lat, plev

def extract_nc_data(files, nc_vari, tmp_array, error_limit):
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
    model_size = 0
    print 'Number of Files', len(files)
    for k in range(len(files)):
        nc = netcdf_file(files[k],'r')
        temp = np.squeeze(nc.variables[nc_variable_name][:])
        temp = temp.copy() # make a copy because temp is a read only version 
        nc.close
        temp[temp>error_limit]= np.nan; # I'd like to make this code genreal for any variable, use dictionary here to assign different tolerances, maybe extract it from the netcdf file?
        temp = np.squeeze(np.nanmean(temp,len(np.shape(temp))-1)); # Take a zonal mean, id like to make this more general as well, automatically pick out the lon dimension?
        print 'temp_file size', np.size(temp,0)
        print 'model_size', model_size
        tmp_array[model_size:model_size+np.size(temp,0),...] = temp
        model_size = model_size + np.size(temp,0)
    return tmp_array
    
def extract_nc_time(files, tmp_array):

    """  Gets time, time from a bunch of netcdf files and concatenates them together

        files: list of pathways to netcdf variables
        tmp_array: An preallocated numpy vector of appropriate length

    Returns:
        A numpy array of all the netcdf files concatenated together

    Raises:
    """

    model_size = 0
    for k in range(len(files)):
        print 'file', k
        nc = netcdf_file(files[k])
        temp = np.squeeze(nc.variables['time'][:])
        temp_copy = temp.copy()
        nc.close
        tmp_array[model_size:model_size+np.size(temp,0)] = temp
        model_size = model_size + np.size(temp,0)
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
        nc = netcdf_file(files[k])
        tmp = np.squeeze(nc.variables[nc_variable_name][:]) # loading dimensions backwards [time, pfull, lat, lon]
        model_size = model_size + np.size(tmp,0);
        nc.close

    return model_size
    
def get_SPOOKIE_filepath(experi,freq,realm,vari,model,mount):
    """Gets a list of the SPOOKIE filepaths from the directory structure I have created,
    
    Directory structure has form SPOOKIE/experi/freq/realm/vari/model/. Operates using a directory mounted up using FUSE OSX
    
    Args:
        experi: experiment (string)
        freq: frequency (string)
        realm: realm (string)
        vari: variable (string)
        model: model (string)
        mount: mountpoint directory name
        
    Returns:
        A list of file pathnames
    """
    location = '/Users/camcairns/' + mount + '/SPOOKIE/'
    files = []
    print "Looking for netcdfs here: \n", location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/r1i1p1/*.nc'
    files = glob.glob(location + str(experi) + '/' + str(freq) + '/' + str(realm) + '/' + str(vari) + '/' + str(model) + '/r1i1p1/*.nc')
    
    return files

def interp_data(lat_old, plev_old, lat_new, plev_new, array):

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

    if plev_old:
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
            if plev_old:
                points =  zip(xold[mask], yold[mask])
                temp = array[m,:,:]
                array_interp[m,:,:] = griddata(points, temp[mask],(xnew,ynew),method='cubic')
            else:
                points = lat_old[mask]
                temp = array[m,:]
                array_interp[m,:] = griddata(points, temp[mask],lat_new,method='cubic')

    return array_interp
        
def write_nc(input_lat, input_latb, input_plev, input_array,  time_array, time_units, time_cal, save_path, nc_vari):

    """ Writes the zonal mean data that has been extracted out as a netcdf file, saves in a new directory structure

        files: list of pathways to netcdf variables
        tmp_array: An preallocated numpy vector of appropriate length

    Returns:
        A numpy array of all the netcdf files concatenated together

    Raises:
    """
    global model, freq, experi, realm, vari

    units_dict = {'ua': 'm/s', 'va': 'm/s', 'uas': 'm/s', 'vas': 'm/s', 'ta': 'K', 'tas': 'K', 'hur': '%', 'hus': '1'}
    # Initiate NETCDF4
    rootgrp = Dataset(save_path, 'w', format='NETCDF4')

    # Set Dimensions
    if input_plev:
        plev = rootgrp.createDimension('plev', 17)
    else:
        pass
    lat = rootgrp.createDimension('lat', 64)
    latb = rootgrp.createDimension('latb', 65)
    time = rootgrp.createDimension('time', None)

    # Set Variables
    times = rootgrp.createVariable('time','f8',('time',))
    if input_plev:
        plev = rootgrp.createVariable('plev','f8',('plev',))
        tmp = rootgrp.createVariable(nc_vari,'f8',('time','plev','lat'))
    else:
        tmp = rootgrp.createVariable(nc_vari,'f8',('time','lat'))
    latitudes = rootgrp.createVariable('latitude','f8',('lat',))
    latitude_bnds = rootgrp.createVariable('latitude_bnds','f8',('latb',))
    # two dimensions unlimited.

    # Attributes for nc file
    import time
    rootgrp.description = 'SPOOKIE data with a zonal mean taken and then interpolated on to a uniform lat-lev grid. This nc_files is for MODEL: ' + model + '  FREQ: ' + freq + '  EXPERIMENT: ' + experi + '  REALM: ' + realm + '  VARI: ' + vari
    rootgrp.history = 'Created ' + time.ctime(time.time()) + '. Original data sourced from DKRZ EUCLIPSE sever.'
    rootgrp.source = 'Cameron Cairns; email: cam.cairns1@gmail.com'
    latitudes.units = 'degrees north'
    if input_plev:
        plev.units = 'hPa'
    else:
        pass
    tmp.units = units_dict[nc_vari]
    times.units = time_units
    times.calendar = time_cal
    
    # Write Data
    latitudes[:] = lat_common
    latitude_bnds[:] = latb_common
    if input_plev:
        plev[:] = plev_common
        tmp[0:model_size,0:len(pfull_common),0:len(lat_common)] = input_array
    else: 
        tmp[0:model_size,0:len(lat_common)] = input_array
    print 'tmp shape after adding data = ', tmp.shape
    
    rootgrp.close()

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
