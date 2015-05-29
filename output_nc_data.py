import numpy as np
import nc_data as nc

def fetch_nc_data(time_length,category,experi, freq, realm, vari, model_list, ensemble=None, mount_dir='mountpoint', start_month=1):
    if time_length in ('max','min'):
        model_size_list = []
        print "The model sizes are:"
        for k, model in enumerate(model_list):
            files = nc.get_filepath(category,experi, freq, realm, vari, model, ensemble=None, mount_dir=mount_dir)
            if files:
                model_size_list.append(nc.find_model_size(files,vari))
            else:
                model_size_list.append(np.nan)
            print '{} : {}'.format(model, model_size_list[k])
        if time_length=='max':
            time_length = max(model_size_list) 
        else:
            time_length = min(model_size_list)
    else:
        pass

    files = nc.get_filepath(category,experi, freq, realm, vari, model_list[0], ensemble, mount_dir=mount_dir)
    output_array, lat, plev, plev_flag = nc.empty_array_generator(files,time_length, model_dim=len(model_list),modulo=12)
    mask_list = []
    for i1, model in enumerate(model_list):
        print model
        files = nc.get_filepath(category,experi, freq, realm, vari, model, ensemble,mount_dir=mount_dir)
        if files:
            model_size = nc.find_model_size(files,vari)
            tmp_array, lat, plev, plev_flag = nc.empty_array_generator(files, model_size)
            tmp_array = nc.extract_nc_data(files, vari, tmp_array, zonal_mean=False)
            time_arg = min(np.size(tmp_array,0),time_length)
            R = nc.round_time(files, model_size,start_month=start_month)
            output_array[:time_arg-R,...,i1] = tmp_array[R:time_arg,...]
        else:
            print "The model %s doesn't exist" % model
            mask_list.append(i1)

    keep_mask = set(range(len(model_list))) - set(mask_list)        
    output_array = output_array[...,list(keep_mask)]
    output_array = nc.reshape_data(output_array,plev_flag)
    output_array = np.squeeze(output_array)
    return output_array
