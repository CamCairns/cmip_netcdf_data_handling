import numpy as np
import nc_data as nc
import os
import glob

def model_intersection(ensemble_list):
    ensemble_models = []
    for dirpath in ensemble_list:
        ensemble_models.append(set(os.listdir(dirpath)))
    for ensemble in ensemble_models:
        ensemble_models[0].intersection_update(ensemble) 
    shared_models = list(ensemble_models[0])
    return shared_models
    
def fetch_nc_data(time_length, experi_list, freq_list, realm_list, vari_list, ensemble=None, mount_dir='mountpoint', start_month=1,verbose=False):
    ensemble_list = nc.generate_ensemble_list(experi_list, freq_list, realm_list, vari_list, mount_dir)
    shared_models = model_intersection(ensemble_list)
    print "Shared model list", shared_models

    if time_length in ('max','min'):
        model_size_list = []
        for experi in experi_list:
            for freq in freq_list:
                for realm in realm_list:
                    for vari in vari_list:
                        for k, model in enumerate(shared_models):
                            files = nc.get_filepath(experi, freq, realm, vari, model, ensemble=None, mount_dir=mount_dir)
                            if files:
                                model_size_list.append(nc.find_model_size(files,vari))
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

    files = nc.get_filepath(experi, freq, realm, vari, shared_models[0], ensemble, mount_dir=mount_dir) # Just getting lat and plev dims (we are assuming all models have shared lat and plev coords
    output_array, lat, plev, plev_flag = nc.empty_array_generator(files,time_length, model_dim=len(shared_models), vari_dim=len(vari_list),modulo=12)

    for experi in experi_list:
        for freq in freq_list:
            for realm in realm_list:
                for i2, vari in enumerate(vari_list):
                    for i1, model in enumerate(shared_models):
                        print experi, vari, model
                        files = nc.get_filepath(experi, freq, realm, vari, model, ensemble,mount_dir=mount_dir)
                        if files:
                            model_size = nc.find_model_size(files,vari)
                            tmp_array, lat, plev, plev_flag = nc.empty_array_generator(files, model_size)
                            tmp_array = nc.extract_nc_data(files, vari, tmp_array, zonal_mean=False)
                            time_arg = min(np.size(tmp_array,0),time_length)
                            R = nc.round_time(files, model_size,start_month=start_month)
                            output_array[:time_arg-R,...,i1,i2] = tmp_array[R:time_arg,...]
                        else:
                            print "The model %s doesn't exist. THAT IS A PROBLEM" % model

#     keep_mask = set(range(len(shared_models))) - set(mask_list)        
#     output_array = output_array[...,list(keep_mask)]
    output_array = nc.reshape_data(output_array,plev_flag)
    output_array = np.squeeze(output_array)
    return output_array
