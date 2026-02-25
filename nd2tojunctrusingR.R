# Testing functionality for JG

# Success looks like: 1. able to source the python script from within R. 2. Able to take nd2 files and convert them to single channel tif MIPs, for quantifying, and also into combined channel MIPs. Combined channels will be uploaded to omero. 3. MIPs look comparable to that created by nikon (sanity check for the MIP processing). 4. Tifs are able to be quantified using the modified vjunctrimport function

# ----------------------------- HOW TO USE THIS CODE -----------------------------

# STEP 1: 
# In the top tab, click code > softwrap long lines. This is important for readability 

# STEP 2: THIS STEP ONLY HAS TO BE DONE ONCE

# You will need to install python onto your computer. And then create a virtual environment that installs packages into that python
# This only has to be done once
install.packages("reticulate")
library(reticulate)

virtualenv_create("r-reticulate")   # create a clean python env
use_virtualenv("r-reticulate", required = TRUE)

py_install(c("numpy", "nd2", "tifffile", "xarray", "pandas"))


# STEP 3:
# You will need a copy of the batchnd2filefunction.py file in your working directory. This is the python script. The package reticulate allows you to source the python script from within Rstudio, without having to use another software

library(reticulate)
use_virtualenv("r-reticulate", required = TRUE) #activate the virtual env

source_python("batchnd2filefunction.py")

# Go to file explorer to the folder where your data is. Copy as path into first argument, changing backslashes to forward slashes
nd2_to_MIPtif("path to folder where original data is",
              "path to a folder you want your two folders to be, and what condition you want the folders to be called")

# for example:
# nd2_to_MIPtif("Y:/Jayden/JI images/Expts for shear paper/5 Dyne 3.0/Channel E&F EVEBcat", 
#              "Y:/Cynthia/Testing python code with JG data/5Dyne3.0ChannelEFEVEBcat", 
 #             overwrite= True)

# This says, in the folder 'Testing python code with JG data, create the folders 5Dyne3.0ChannelEFEVEBcat_for OMERO, and 5Dyne3.0ChannelEFEVEBcat_singlechannelTIF