import numpy as np
import nc_data as ncd
import os
import glob

def generate_shared_models(experi_list, freq_list, realm_list, vari_list, root_list = None):
    ensemble_list = ncd.generate_ensemble_list(experi_list, freq_list, realm_list, vari_list)
    shared_models = ncd.model_intersection(ensemble_list)
    if root_list:
        shared_models = list(set(shared_models).intersection(set(root_list)))
    return shared_models

def fetch_nc_data(time_length, experi_list, freq_list, realm_list, vari_list, root_list = None, start_month=1, verbose=False):
    """Fetches netcdf data from SPOOKIE/AMIP interpolated directory structure.

 Outputs this data as a numpy array with form [month year plev lat model variable]. 
 
 The function does a couple of things  
    1) The time length of the data array is specified by the user (see args: time_length)
    2) The models in the output array are the common elements of the set of models shared between experiment/variable directories (as well as a root list if it is specified)
    3) Detects the start month of a model and anchors each model to a shared start month (default = January)
    4) pads out the array with NaNs so that it can be reshaped into [month year] dims

All combinations of the various input lists are formed 
For example:
        experi_list = ['SPOOKIE_interp/convoffamip','AMIP_interp/amip']
        freq_list = ['mon']
        realm_list = ['atmos']
        vari_list = ['ta','va'] 
    
would find the 4 combinationsL
        'SPOOKIE_interp/convoffamip/mon/atmos/ta'
        'SPOOKIE_interp/convoffamip/mon/atmos/ua'
        'AMIP_interp/amip/mon/atmos/ta'
        'AMIP_interp/amip/mon/atmos/ua'
        
    Args:
        time_length: integer, 'max' or 'min'. Sets the time_length of the eventual output_array
             'max' ==> an array size determined by the maximum time length of shared models in the available ensemble (padded up to the next multiple of 12)
             'min' ==> an array size determined by the minimum time length of shared models in the available ensemble (padded up to the next multiple of 12)
             <integer> ==> an array size with time_length <integer> (padded up to the next multiple of 12)
        experi_list: list of experiments for example ['SPOOKIE/convoffamip']
        freq_list: list of freq terms
        realm_list: list of realm terms
        vari_list: list of vari terms
        root_list: a user specified list of models to intersect with the available ensemble list. If nothing is specified for root_list the full available ensemble will be used
        start_month (optional): The anchor month from which each time series begins from (default = 1 = jan)
        verbose (optional): if True, ouputs some helpful print statements

    Returns:
        output_array: A numpy array with dimensions [month, year, plev, lat, common_models, vari]
        shared_models: A list of the shared models comprising output_array. The list is in the order they appear in the model dimension of output_array
    """
    shared_models = generate_shared_models(experi_list, freq_list, realm_list, vari_list, root_list = root_list)
    print "Shared model list", shared_models

    time_length = ncd.get_time_dim(experi_list, freq_list, realm_list, vari_list, shared_models, time_length=time_length)
    files = ncd.get_filepath(experi_list[0], freq_list[0], realm_list[0], vari_list[0], shared_models[0]) # Just getting lat and plev dims (we are assuming all models have shared lat and plev coords
    plev, lat, lon, plev_flag, latb, lonb = ncd.load_coord_data(files)
    time_length = ncd.modulo_padding(time_length,12)
    output_array = ncd.empty_array_generator([time_length, len(plev), len(lat), len(lon), len(shared_models), len(experi_list), len(vari_list)])
    for i5, experi in enumerate(experi_list):
        for freq in freq_list:
            for realm in realm_list:
                for i2, vari in enumerate(vari_list):
                    for i1, model in enumerate(shared_models):
                        files = ncd.get_filepath(experi, freq, realm, vari, model, verbose=verbose)
                        if files:
                            model_size = ncd.find_model_size(files,vari)
                            plev, lat, lon, plev_flag, latb, lonb = ncd.load_coord_data(files)
                            tmp_array = ncd.empty_array_generator([model_size, len(plev), len(lat), len(lon)])
                            tmp_array = ncd.extract_nc_data(files, vari, tmp_array, zonal_mean=False)
                            time_arg = min(np.size(tmp_array,0),time_length)
                            R = ncd.round_time(files, model_size,start_month=start_month)
                            output_array[:time_arg-R,...,i1,i5,i2] = tmp_array[R:time_arg,...]
                        else:
                            print "The model %s doesn't exist. THAT IS A PROBLEM" % model
    output_array = ncd.reshape_data(output_array,plev_flag)
    print 'output array shape', np.shape(output_array)
    return output_array, lat, plev, latb, lonb, shared_models
