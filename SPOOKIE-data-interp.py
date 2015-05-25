# Main Program

model_list = ['CanAM4','CNRM-AM6', 'CNRM-CM5', 'MIROC5', 'HadGEM2-A', 'MPI-ESM-LR', 'MRI-CGCM3', 'GFDL-HIRAM', 'CESM-CAM5.1-FV2']
freq_list = ['mon']
experi_list = ['convoffamip4K'] #,'convoffamip4xCO2'] #['convoffamip'] #,
realm_list = ['atmos']
vari_list = ['ua','va','ta','hur','hus']
mount = 'mountpoint3'

for i1, experi in enumerate(experi_list):
    for i2, freq in enumerate(freq_list):
        for i3, realm in enumerate(realm_list):
            for i4, vari in enumerate(vari_list):
                for i5, model in enumerate(model_list):
                    print experi, freq, realm, vari, model
                    files = get_SPOOKIE_filepath(experi, freq, realm, vari, model, mount)
                    if files:
                        files = get_SPOOKIE_filepath(experi, freq, realm, vari, model, mount)
                        plev_common = [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000]
                        latb_common= np.linspace(-90,90,65);
                        lat_common = 0.5*(3*latb_common[0]-latb_common[1])+np.cumsum(np.diff(latb_common));

                        time_init= np.empty([1]) # initalise newtemp with the correct dimensions
                        time_units = []
                        time_cal = []
                    
                        lat, plev = load_coord_data(files)
                    
                        model_size = find_model_size(files, vari) # find the total time length of the model, use to preallocate a numpy array
                        print 'model_size', model_size
                    
                        if plev:
                            tmp_array = np.empty([model_size,len(plev),len(lat)])*np.nan;
                        else:
                            plev_common = None
                            tmp_array = np.empty([model_size,len(lat)])*np.nan;
                    
                        #FIELD DATA EXTRACTION AND CONCATENATION    
                        tmp_array = extract_nc_data(files, vari, tmp_array, model_size, 1.0e8);

                        #INTERPOLATE ONTO COMMON GRID
                        tmp_array_interp = interp_data(lat, plev, lat_common, plev_common, tmp_array)                    
                    
                        # TIME DATA EXTRACTION
                        tmp_array = np.empty([model_size])*np.nan;
                        time_array = extract_nc_time(files, tmp_array)
                        nc = netcdf_file(files[0])
                        time_units = np.append(time_units,nc.variables['time'].units)
                        time_cal = np.append(time_cal,nc.variables['time'].calendar)
                        nc.close
                    
                        # Make directory path                    
                        save_path = '/Users/camcairns/' + mount + '/SPOOKIE_interp/' + experi + '/' + freq + '/' + realm + '/' + vari + '/' + model
                        mkdir_p(save_path)
                        # Write nc files
                        nc_file = vari + '_' + freq + '_' + experi + '_SPOOKIE_interp.nc'
                        write_nc(lat_common, latb_common, plev_common, tmp_array_interp, time_array, time_units, time_cal, save_path + '/' + nc_file, vari)
                    else: 
                        pass
