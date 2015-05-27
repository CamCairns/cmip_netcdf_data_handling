import numpy as np
import nc_data as nc

category = 'SPOOKIE_interp'
experi = 'convoffamip'
freq = 'mon'
realm = 'atmos'
vari =  'ua'
model_list = ['CanAM4','CNRM-AM6', 'CNRM-CM5','NotInThere']#, 'MIROC5', 'HadGEM2-A', 'MPI-ESM-LR', 'MRI-CGCM3', 'GFDL-HIRAM', 'CESM1-CAM5']
mount = 'mountpoint3'

def fetch_nc_data(time_length,category,experi, freq, realm, vari, model_list, ensemble=None, mount_dir=mount):
    if time_length in ('max','min'):
        model_size_list = []
        print "The model sizes are:"
        for k, model in enumerate(model_list):
            files = nc.get_filepath(category,experi, freq, realm, vari, model, ensemble=None, mount_dir=mount)
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

    files = nc.get_filepath(category,experi, freq, realm, vari, model_list[0], ensemble=None, mount_dir=mount)
    output_array, lat, plev, plev_flag = nc.empty_array_generator(files,time_length, model_dim=len(model_list))
    mask_list = []
    for i1, model in enumerate(model_list):
#        print "Loading %s/%s/%s/%s/%s/%s" % (category, experi, freq, realm, vari, model)
        files = nc.get_filepath(category,experi, freq, realm, vari, model, ensemble=None,mount_dir=mount)
        if files:
            model_size = nc.find_model_size(files,vari)
            tmp_array, lat, plev, plev_flag = nc.empty_array_generator(files, model_size)
            tmp_array = nc.extract_nc_data(files, vari, tmp_array, zonal_mean=False)
            time_arg = min(np.size(tmp_array,0),time_length)
            output_array[:time_arg,...,i1] = tmp_array[:time_arg,...]
        else:
            print "The model %s doesn't exist" % model
            mask_list.append(i1)

    keep_mask = set(range(len(model_list))) - set(mask_list)        
    output_array = output_array[...,list(keep_mask)]
    return output_array