# Main Program
from nc_data import *

model_list = ['CanAM4','CNRM-AM6', 'CNRM-CM5', 'MIROC5', 'HadGEM2-A', 'MPI-ESM-LR', 'MRI-CGCM3', 'GFDL-HIRAM', 'CESM1-CAM5','CESM-CAM5.1-FV2']
freq_list = ['mon']
experi_list = ['amip4K','amip4xCO2', 'aqua4xCO2','aqua4K', 'aquaControl']#,'convoffaqua4K','convoffaquaControl'] #['amip4K', 'amip4xCO2', 'aqua4xCO2', 'aqua4K', 'aquaControl'] #[,'convoffamip4xCO2']
exp_dir = 'AMIP'
realm_list = ['atmos']
vari_list = ['uas','vas','tas']#,'vas','tas']#['ua','uas','va','vas','ta','tas','hur','hus','zg']
ensemble = 'r1i1p1'

for i1, experi in enumerate(experi_list):
    for i2, freq in enumerate(freq_list):
        for i3, realm in enumerate(realm_list):
            for i4, vari in enumerate(vari_list):
                for i5, model in enumerate(model_list):
                    print experi, freq, realm, vari, model
                    files = get_filepath(exp_dir + '/' + experi, freq, realm, vari, model)
                    if files:
                        if not os.path.isdir(os.path.join(location, exp_dir + '_interp', experi, freq, realm, vari, model)):
                            plev_common = [100000, 92500, 85000, 70000, 60000, 50000, 40000, 30000, 25000, 20000, 15000, 10000, 7000, 5000, 3000, 2000, 1000]
                            latb_common= np.linspace(-90,90,65);
                            lat_common = 0.5*(3*latb_common[0]-latb_common[1])+np.cumsum(np.diff(latb_common));

                            time_init= np.empty([1]) # initalise newtemp with the correct dimensions
                            time_units = []
                            time_cal = []
                    
                            model_size = find_model_size(files, vari) # find the total time length of the model, use to preallocate a numpy array
                        
                            # PREALLOCATE EMPTY ARRAY
                            plev, lat, lon, plev_flag, latb, lonb = load_coord_data(files)
                            tmp_array = empty_array_generator([model_size, len(plev), len(lat)])

                            #FIELD DATA EXTRACTION AND CONCATENATION    
                            tmp_array = extract_nc_data(files, vari, tmp_array);

                            #INTERPOLATE ONTO COMMON GRID
                            tmp_array_interp = interp_data(plev,plev_flag,lat,plev_common,lat_common,tmp_array)
                    
                            # TIME DATA EXTRACTION
                            time_vector, time_units, time_cal = extract_nc_time(files, model_size)

                            # Make directory path
                            save_path = os.path.join(location, exp_dir + '_interp', experi, freq, realm, vari, model, ensemble)
                            mkdir_p(save_path)
                            # Write nc files
                            nc_file = vari + '_' + freq + '_' + model.replace('.','-') + '_' + experi + '_' + exp_dir + '_interp.nc'
                            write_nc(lat_common, latb_common, plev_common, plev_flag, tmp_array_interp, time_vector, time_units, time_cal, save_path + '/' + nc_file, model_size, experi, freq, realm, vari, model)
                            print 'Netcdf file written for /%s/%s/%s/%s/%s' % (experi, freq, realm, vari, model)
                        else:
                            print 'That path already exists'
                    else: 
                        print 'No netcdf files found in that location!'
