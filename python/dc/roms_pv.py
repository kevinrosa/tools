def roms_pv(fname,tindices):

    import h5py
    u = dc_roms_read_data(fname,'u')
    v = dc_roms_read_data(fname,'v')
    
