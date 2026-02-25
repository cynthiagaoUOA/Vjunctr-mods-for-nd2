library(tools)
library(tidyverse)


# Import function -------------------------------------------------------


# For Jayden's data, where file names are 01_A01_ch2.tif, every loaded in folder will be it's own condition (equivalent to my well),
# and each "well" or A1 etc will be the "fields". Keep the wording of fields I think. 

JGimport_maxIPtif = function(file_path, lookup = NULL){
  
  # Grab the file list
  filelist = list.files(file_path, full.names = TRUE, pattern = ".tif") #grabs a list of all filenames in the provided tiff folder
  file_directory = data.frame(path = filelist)
  
  # Grab the metadata from the file and split it out
  file_directory$file = basename(file_path_sans_ext(file_directory$path)) # file column that ignore full path just pastes the file name, eg 02_B02_0001_ch0.tif

  file_directory$field <- match(substr(file_directory$file, 4, 4), LETTERS) #Takes the letter, turns it into number, as field
  file_directory$shearchannel <- substr(file_directory$file, 5, 5)
  
  file_directory$channeltemporary = substr(file_directory$file,9,9) %>% as.numeric() # channel as in colours
  file_directory$channeltemporary = file_directory$channeltemporary +1 #so field starts from 1
  file_directory$channel = paste("CH_",file_directory$channeltemporary, sep = "")
  
  file_directory$file = NULL #don't need this column anymore, all key info split into other columns
  file_directory$channeltemporary = NULL
  
  
  file_directory = as.tibble(file_directory)
  
  if(is.null(lookup))
  {
    return(file_directory)
  }
  
  file_directory = file_directory %>% pivot_wider(names_from = channel, values_from = path)
  
  # Attach the lockup table data
  file_directory = right_join(file_directory, lookup, by = "shearchannel", relationship="many-to-many")
  
  
  # Move the colnames where they need to be
  file_directory = file_directory %>%
    mutate(ch_dapi = case_when(ch_dapi == 1 ~ CH_1,
                               ch_dapi == 2 ~ CH_2,
                               ch_dapi == 3 ~ CH_3,
                               ch_dapi == 4 ~ CH_4)) %>%
    mutate(ch_actin = case_when(ch_actin == 1 ~ CH_1,
                                ch_actin == 2 ~ CH_2,
                                ch_actin == 3 ~ CH_3,
                                ch_actin == 4 ~ CH_4)) %>%
    mutate(ch_antibody = case_when(ch_antibody == 1 ~ CH_1,
                                   ch_antibody == 2 ~ CH_2,
                                   ch_antibody == 3 ~ CH_3,
                                   ch_antibody == 4 ~ CH_4)) %>%
    select(-CH_1, -CH_2, -CH_3, -CH_4)
  
  
 # Want to match conditions, potentially up to six shear channels
  
  return(file_directory)
  
}



# EXAMPLE TEMPLATE --------------------------------------------------------

## COde tested w JG data, runs. Copied below as a template

library(tidyverse)
library(BiocManager)
library(EBImage)
library(data.table)
library("shiny")
library("bslib")
library(progressr)
library(doFuture)
library(patchwork)
library(gglm)

# load the enviornment that contains all installed python packages
library(reticulate)
use_virtualenv("r-reticulate", required = TRUE) #activate the virtual env

source_python("nd2tojunctr2.0.py") #allows nd2_to_MIPtif("Y:/Jayden/JI images/Expts for shear paper/5 Dyne 3.0/Channel E&F CE_VEBcat",

source("JI_nikon_importfunction.R") #source function for importing tifs into labeled format
source("vjunctur_functions_v1.R") # the rest of JH's vjunctr code and functions


#####

# Turn nd2 files into tifs
nd2_to_MIPtif("Y:/Jayden/JI images/Expts for shear paper/5 Dyne 3.0/Channel E&F CE_VEBcat",
              "Y:/Jayden/JI images/Expts for shear paper/5 Dyne 3.0/Channel E&F", overwrite= FALSE)


##### import data and add plate map

# 488 VE-cad goat, Dapi, B-cat then actin

platemap = tribble(~shearchannel, ~ch_dapi, ~ch_actin, ~ch_antibody, ~name_antibody, ~sample,
                   "1", 2, 4, 1, "ve-cadherin", "shear no melanoma", #match shearchannel location to sample/condition
                   "1", 2, 4, 3, "b-catenin", "shear no melanoma", #each line a different antibody. Dapi and actin in all
                   
                   "2", 2, 4, 1, "ve-cadherin", "static no melanoma",
                   "2", 2, 4, 3, "b-catenin", "static no melanoma")

data<- JGimport_maxIPtif("Y:/Jayden/JI images/Expts for shear paper/5 Dyne 3.0/Channel E&F_singlechannelTIF", platemap)
data$rep = "run1"       #creates a column tagging this dataset as replicate 1. This way when combining runs later, can differentiate

# Quantify eg VE-cad
vecad <- data %>% filter(name_antibody== "ve-cadherin")

segment_and_quant_i(vecad)

