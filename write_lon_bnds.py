import numpy as np
from netCDF4 import Dataset
import nc_data as ncd

def write_field_value(file_path,string_name,input_value):
    nc = Dataset(file_path,'r+')
    nc.variables[string_name][:,:]=input_value
    nc.close

# if __name__ == '__main__':
def write_lon_bnds():
'''
Creates a variable lon_bnds and populates it with data from the corresponding AMIP experiment. Called from write_lonb_latb.py. 

NOTE: write_lon_bnds.py must be called BEFORE write_lat_bnds.py as this function creates a bnds variable that is needed here 

'''
    model_list = ['MPI-ESM-LR'] 
    freq_list = ['mon']
    experi_list = ['SPOOKIE/convoffamip']#,'SPOOKIE/convoffamip4xCO2']#,'SPOOKIE/convoffaqua4K','SPOOKIE/convoffaqua4xCO2','SPOOKIE/convoffaquaControl'] 
    AMIP_experi_list = ['AMIP/amip']#, 'AMIP/amip4xCO2']#, 'AMIP/aqua4xCO2', 'AMIP/aqua4K', 'AMIP/aqua4xCO2', 'AMIP/aquaControl']
    realm_list = ['atmos']
    vari_list = ['uas','va','vas','ta','tas','hur','hus']

    for i1, experi in enumerate(experi_list):
        for freq in freq_list:
            for realm in realm_list:
                for vari in vari_list:
                    files_AMIP = ncd.get_filepath(AMIP_experi_list[i1], freq, realm, vari, 'MPI-ESM-LR')
                    print files_AMIP[0]
                    nc = Dataset(files_AMIP[0],'r')
                    input_lonb = nc.variables['lon_bnds'][:]
                    nc.close
                    files_SPOOKIE = ncd.get_filepath(experi, freq, realm, vari, 'MPI-ESM-LR')
                    if files_SPOOKIE:
                        print files_SPOOKIE[0]
                        nc = Dataset(files_SPOOKIE[0],'r+')
                        # Set Dimensions
                        bnds = nc.createDimension('bnds', np.size(input_lonb,1))
                        # Set Variables
                        longitude_bnds = nc.createVariable('lon_bnds','f8',('lon','bnds'))
                        nc.close
                        write_field_value(files_SPOOKIE[0],'lon_bnds',input_lonb)
                    else: 
                        print 'That file does not exist!'