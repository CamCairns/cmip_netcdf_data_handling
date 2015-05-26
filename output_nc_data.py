from nc_data import  nc_data as cn
time_length =  # max_length, min_length or some integer

for i1, experi in enumerate(experi_list):
    for i2, freq in enumerate(freq_list):
        for i3, realm in enumerate(realm_list):
            for i4, vari in enumerate(vari_list):
                for i5, model in enumerate(model_list):
                    print experi, freq, realm, vari, model
                    files = get_SPOOKIE_filepath(category,experi, freq, realm, vari, model, mount)
                    if files:
                        
                        #FIELD DATA EXTRACTION AND CONCATENATION    
                        tmp_array = extract_nc_data(files, nc_vari, tmp_array, zonal_mean=False)