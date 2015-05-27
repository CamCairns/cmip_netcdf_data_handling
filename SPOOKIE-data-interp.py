# Main Program
from nc_data import *

model_list = ['CanAM4','CNRM-AM6', 'CNRM-CM5', 'MIROC5', 'HadGEM2-A', 'MPI-ESM-LR', 'MRI-CGCM3', 'GFDL-HIRAM', 'CESM1-CAM5'] #'CESM-CAM5.1-FV2'
freq_list = ['mon']
experi_list = ['aqua4xCO2', 'aquaControl'] #['convoffaqua4K','convoffaqua4xCO2','convoffaquaControl'] #,'convoffamip4xCO2'] #
realm_list = ['atmos']
vari_list = ['ua','uas','va','vas','ta','tas','hur','hus']#
category = 'AMIP' # SPOOKIE
mount = 'mountpoint3'

for i1, experi in enumerate(experi_list):
    for i2, freq in enumerate(freq_list):
        for i3, realm in enumerate(realm_list):
            for i4, vari in enumerate(vari_list):
                for i5, model in enumerate(model_list):
                    print experi, freq, realm, vari, model
                    files = get_filepath(category,experi, freq, realm, vari, model, mount_dir = mount)
                    if files:
                        plev_common = [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000]
                        latb_common= np.linspace(-90,90,65);
                        lat_common = 0.5*(3*latb_common[0]-latb_common[1])+np.cumsum(np.diff(latb_common));

                        time_init= np.empty([1]) # initalise newtemp with the correct dimensions
                        time_units = []
                        time_cal = []
                    
                        model_size = find_model_size(files, vari) # find the total time length of the model, use to preallocate a numpy array
                        
                        # PREALLOCATE EMPTY ARRAY
                        tmp_array, lat, plev, plev_flag = empty_array_generator(files, model_size)

                        #FIELD DATA EXTRACTION AND CONCATENATION    
                        tmp_array = extract_nc_data(files, vari, tmp_array);

                        #INTERPOLATE ONTO COMMON GRID
                        tmp_array_interp = interp_data(lat, plev, plev_flag, lat_common, plev_common, tmp_array)                    
                    
                        # TIME DATA EXTRACTION
                        tmp_array = np.empty([model_size])*np.nan;
                        time_array = extract_nc_time(files, tmp_array)
                        nc = Dataset(files[0])
                        time_units = np.append(time_units,nc.variables['time'].units)
                        time_cal = np.append(time_cal,nc.variables['time'].calendar)
                        nc.close
                    
                        # Make directory path                    
                        save_path = '/Users/camcairns/' + mount + '/' + category + '_interp/' + experi + '/' + freq + '/' + realm + '/' + vari + '/' + model
                        mkdir_p(save_path)
                        # Write nc files
                        nc_file = vari + '_' + freq + '_' + experi + '_' + category + '_interp.nc'

                        write_nc(lat_common, latb_common, plev_common, plev_flag, tmp_array_interp, time_array, time_units, time_cal, save_path + '/' + nc_file, model_size, experi, freq, realm, vari, model)
                        print 'Netcdf file written for /%s/%s/%s/%s/%s' % (experi, freq, realm, vari, model)
                    else: 
                        print 'No netcdf files found in that location!'
