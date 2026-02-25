import pandas as pd
import numpy as np

import nd2 
import tifffile
from tifffile import imwrite
import xarray
from pathlib import Path




def CGnd2tojunctr(input_folder, output_folder):  # defines a function where I can specifiy a folder of nd2 files and it will turn them into single channel MIP TIFs

    input_folder = Path(input_folder) # define what folder of nd2 files to read in
    output_folder = Path(output_folder)
    output_folder.mkdir(exist_ok=True) # if folder doesn't exist, create it. If it does, this code prevents script from breaking

    for nd2_file in input_folder.glob("*.nd2"):  # for every ND2 file in the input folder, loop
        print(f"Processing {nd2_file.name}")    #progress message

        with nd2.ND2File(nd2_file) as f:         # make the file an object called f
            x=f.to_xarray()                     # turn f into an array
   #Maximum intensities
            if 'Z' in x.dims:
                mip = x.max(dim='Z')            # if stacks, make MIP
            else:
                raise ValueError(f"{nd2_file.name} has no Z axis — cannot compute MIP.") #otherwise print error msg
    
    # Split channels 
            if 'C' in mip.dims:
                for c in range(mip.sizes['C']):     # For every channel
                    img = np.asarray(mip.isel(C=c)) # Turn into array
                    out_path = output_folder / f"{nd2_file.stem}_ch{c}.tif" # specify a location and file name. Using original base name, with channel number added
                    imwrite(out_path, img)          # at this output path, create tif image
            else:
                img = np.asarray(mip)
                out_path = output_folder / f"{nd2_file.stem}_MIP.tif" # if no channels just make original into tif
                imwrite(out_path, img)
    print(f"Finished {nd2_file.name}\n")


    
def nd2_to_MIPtif(input_folder, output_folder, overwrite= False):  # defines a function where I can specifiy a folder of nd2 files and it will turn them into single channel MIP TIFs

    input_folder = Path(input_folder) # define what folder of nd2 files to read in
    output_folder = Path(output_folder)

 # create two folders using the output folder specified. 

    parent= output_folder.parent # I want to create two output folders, so this is specifying the parent, in order to define the subfolders
    base = output_folder.name # defining the base name as whatever I put into the function

    omero_dir = parent / f"{base}_MIPs_forOMERO" # create first subfolder. 
    omero_dir.mkdir(parents=True, exist_ok=True) # if folder doesn't exist, create it. If it does, this code prevents script from breaking
    single_dir = parent / f"{base}_singlechannelTIF" # same as above
    single_dir.mkdir(parents=True, exist_ok=True)

 
    for nd2_file in input_folder.glob("*.nd2"):  # for every ND2 file in the input folder, do the following
        print(f"Processing {nd2_file.name}")    #progress message

        with nd2.ND2File(nd2_file) as f:         # make the file an object called f
            x=f.to_xarray()                     # turn into an array
   #Make everything into MIPS
            if 'Z' in x.dims:
                mip = x.max(dim='Z')            # if stacks, make MIP
            else:
                raise ValueError(f"{nd2_file.name} has no Z axis — cannot compute MIP.") #otherwise print error msg
    
    # # OUTPUT 1: In the singlechannelTIF folder make single channel TIFs   
            if 'C' in mip.dims: #given that there are multiple channels in the file,
                for c in range(mip.sizes['C']):     # For every channel
                    img = np.asarray(mip.isel(C=c)) # Turn into array for individual channels
                    out_path = single_dir / f"{nd2_file.stem}_ch{c}.tif" #specifies outputs as tifs with corresponding channel number in the TIF folder
                    if overwrite or not out_path.exists(): #If either we specified to overwirte, or if file doesn't exist yet, then
                        imwrite(out_path, img) # create that tif

            else:
                img = np.asarray(mip) # If there are not multiple channels then just leave as a MIP 
                out_path = single_dir / f"{nd2_file.stem}.tif"
                if overwrite or not out_path.exists():
                    imwrite(out_path, img)

    # OUTPUT 2: Combined channel tifs
            if "C" in mip.dims: # If there are multiple channels in the file
                arr_yxc = np.asarray(mip.transpose("Y", "X", "C")) #create array, grayscale Tif OME-TIFF. 
                out_ome = omero_dir / f"{nd2_file.stem}.ome.tif" # Put these outputs in the omero folder
                if overwrite or not out_ome.exists():
                    imwrite(out_ome, arr_yxc, metadata={"axes": "YXC"}, photometric="minisblack") # preserves metadata for omero/image readers apparently
            else:
                arr_yx = np.asarray(mip)
                out_gray = omero_dir / f"{nd2_file.stem}.tif" # Single channel image just create the Tif
                if overwrite or not out_gray.exists():
                    imwrite(out_gray, arr_yx)
        print(f"Finished {nd2_file.name}\n")

