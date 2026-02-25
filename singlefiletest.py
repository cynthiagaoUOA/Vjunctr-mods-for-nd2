print("hello")
import pandas as pd
import numpy as np

import nd2 
import tifffile
from tifffile import imwrite
import xarray


with nd2.ND2File("01_B02_0000.nd2") as f:
    print(f.sizes)
    print(f.dtype)
    
    x = f.to_xarray()

    # ensure Z exists
    if 'Z' in x.dims:
        mip = x.max(dim='Z')  # now dims likely include C, Y, X
    else:
        raise ValueError("This ND2 has no Z axis (no stack), so MIP doesn't apply.")


    # If multiple channels, save one per channel
    if 'C' in mip.dims:
        for c in range(mip.sizes['C']):
            img = np.asarray(mip.isel(C=c))          # convert to a plain array
            imwrite(f"example_MIP_ch{c}.tif", img)   # write TIFF
    else:
        img = np.asarray(mip)
        imwrite("example_MIP.tif", img)



print ("hi")