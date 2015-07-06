import nc_data as ncd

'''
Function creates lonb and latb variables and populates the variable with data drom the corresponding AMIP experiment for the MPI-ESM-LR SPOOKIE
experiment.
 
These are INCREDIBLY delicate functions. Sometimes they write, sometimes they don't. They seem to work in there current configuration 
but when I tried to move the coded into one .py file they broke again in a strange way (only half the lonb values were written).

I really don't understand why they are so finicky and tbh I really don't care. They do the (one off) job.
'''
if __name__ == '__main__':
# The order of these two functions is important
    ncd.write_lon_bnds()
    ncd.write_lat_bnds() 