from netCDF4 import Dataset
import nc_data as ncd

model_list = ['MPI-ESM-LR'] #['CanAM4','CNRM-AM6', 'CNRM-CM5', 'MIROC5', 'HadGEM2-A', 'MPI-ESM-LR', 'MRI-CGCM3', 'GFDL-HIRAM', 'CESM1-CAM5','CESM-CAM5.1-FV2']
freq_list = ['mon']
experi_list = ['SPOOKIE/convoffamip4K','SPOOKIE/convoffamip4xCO2']#,'SPOOKIE/convoffaqua4K','SPOOKIE/convoffaqua4xCO2','SPOOKIE/convoffaquaControl'] 
AMIP_experi_list = ['AMIP/amip4K', 'AMIP/amip4xCO2']#, 'AMIP/aqua4xCO2', 'AMIP/aqua4K', 'AMIP/aqua4xCO2', 'AMIP/aquaControl']
exp_dir = 'SPOOKIE'
realm_list = ['atmos']
vari_list = ['ua','uas','va','vas','ta','tas','hur','hus']
mount = 'mountpoint3'

for i1, experi in enumerate(experi_list):
    for freq in freq_list:
        for realm in realm_list:
            for vari in vari_list:
                files_AMIP = ncd.get_filepath(AMIP_experi_list[i1], freq, realm, vari, 'MPI-ESM-LR')
                print files_AMIP[0]
                nc = Dataset(files_AMIP[0],'r')
#                 time_vect = nc.variables['time'][:]
#                 time_units = nc.variables['time'].units
                time_vect = nc.variables['time'][:]
                time_units = nc.variables['time'].units
                nc.close
                files_SPOOKIE = ncd.get_filepath(experi, freq, realm, vari, 'MPI-ESM-LR')
                if files_SPOOKIE:
                    print files_SPOOKIE[0]
                    nc = Dataset(files_SPOOKIE[0],'r+')
                    nc.variables['time'][:]=time_vect
                    nc.variables['time'].units=time_units
                    nc.close
                else: 
                    print 'That file does not exist!'