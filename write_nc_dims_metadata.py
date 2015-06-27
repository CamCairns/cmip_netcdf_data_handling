from netCDF4 import Dataset
import nc_data as ncd

model_list = ['MPI-ESM-LR'] 
freq_list = ['mon']
experi_list = ['SPOOKIE/convoffamip4K','SPOOKIE/convoffamip4xCO2']#,'SPOOKIE/convoffaqua4K','SPOOKIE/convoffaqua4xCO2','SPOOKIE/convoffaquaControl'] 
AMIP_experi_list = ['AMIP/amip4K', 'AMIP/amip4xCO2']#, 'AMIP/aqua4xCO2', 'AMIP/aqua4K', 'AMIP/aqua4xCO2', 'AMIP/aquaControl']
exp_dir = 'SPOOKIE'
realm_list = ['atmos']
vari_list = ['ua','uas','va','vas','ta','tas','hur','hus']

for i1, experi in enumerate(experi_list):
    for freq in freq_list:
        for realm in realm_list:
            for vari in vari_list:
                files_AMIP = ncd.get_filepath(AMIP_experi_list[i1], freq, realm, vari, 'MPI-ESM-LR')
                print files_AMIP[0]
                nc = Dataset(files_AMIP[0],'r')
                input_latb = nc.variables['lat_bnds'][:]
                input_lonb = nc.variables['lon_bnds'][:]
                nc.close
                files_SPOOKIE = ncd.get_filepath(experi, freq, realm, vari, 'MPI-ESM-LR')
                if files_SPOOKIE:
                    print files_SPOOKIE[0]
                    nc = Dataset(files_SPOOKIE[0],'r+')
                    # Set Dimensions
                    latb = rootgrp.createDimension('lat_bnds', len(input_latb))
                    lonb = rootgrp.createDimension('lon_bnds', len(input_lonb))
                    # Set Variables
                    latitude_bnds = rootgrp.createVariable('lat_bnds','f8',('latb',))
                    longitudes_bnds = rootgrp.createVariable('lon_bnds','f8',('lonb',))
                    # Write Data
                    latitude_bnds[:] = input_latb
                    longitudes_bnds[:] = input_lonb
                    nc.close
                else: 
                    print 'That file does not exist!'