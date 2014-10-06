def dc_roms_read_data(folder, varname, tindices = [1 Inf], volume = [], stride=[1 1 1 1], grd):
    import numpy as np
    import h5py
    import glob
    import os
    import pyroms

    if os.path.isdir(folder):
        files = glob.glob(folder + '/*avg*');
    
    for ii in range(len(files)):
        if os.path.isdir(folder):
            file = h5py.File(files[ii],'r')
        else:
            file = h5py.File(folder,'r')
            
        var = file[varname]
        nt = var.shape[0]

        if ii == 0:
            ndim = len(var.shape)
            
            if ndim ~= 1:
                # get grid stuff here
                
            else:
                vol = []
                xax = []
                yax = []
                zax = []
                stride[1:] = []
                
            # allocate data here
            size = np.zeros(ndim,np.int32)
            # this must be changed based on vol
            size[-1] = var.shape[-1]
            size[-2] = var.shape[-2]
            size[0] = tindices[1]
            varout = np.zeros(size)

        # process tindices here
        
        if np.isinf(tindices[1]):
            tnew[1] = Inf

        stride[-1] = 1
        
        # Case 1: if requested data is not in this file, skip to next
        if tnew[0] > nt:
            tindices = tindices - nt
            continue
        else:
            # Case 2 : requested data spills into next file
            if tnew[0] <= nt && tnew[1] > nt:
                # change ending for current read
                tnew[1] = nt
                # set next read to start from beginning of new file
                tindices[0] = 1
                tindices[1] -= nt
            else:
                # Case 3 : requested data finishes in current file
                if tnew[1] <= nt && ~np.isinf(tindices[1]):
                    quitflag = 1
                    
        # save data
        if ndim == 4:
            varout[k:,:,:,:] = var[XXX,:,:,:]
        else:
            if ndim == 3:
            varout[k:,:,:] = var[XXX,:,:]
        
        else:
            if ndim == 2:
            varout[k:,:] = var[XXX,:,:]
        else:
            varout[k:] = var[XXX]

    
    return varout


                
