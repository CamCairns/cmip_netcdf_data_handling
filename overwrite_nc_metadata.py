from netCDF4 import Dataset

files = '/Volumes/MyBook/SPOOKIE/convoffamip/mon/atmos/tas/MPI-ESM-LR/r1i1p1/tas_Amon_MPI-ESM-LR_convoffamip_r1i1p1_19790101-20081231.nc'
nc = Dataset(files,'r+')
correct_time = 'days since 1979-1-1 00:00:00'
nc.variables['time'].units=correct_time
nc.close