library(tools)
library(tidyverse)

# Cynthia - wells ---------------------------------------------------------

#prereq - input should be a folder where each tif is an individual channel for one field within a well.

CGimport_maxIPtif = function(file_path, lookup = NULL){
  
  # Grab the file list
  filelist = list.files(file_path, full.names = TRUE, pattern = ".tif") #grabs a list of all filenames in the provided tiff folder
  file_directory = data.frame(path = filelist)
  
  # Grab the metadata from the file and split it out
  file_directory$file = basename(file_path_sans_ext(file_directory$path)) # file column that ignore full path just pastes the file name, 
  #eg 02_B02_0001_ch0.tif
  
  file_directory$well = substr(file_directory$file,4,6) #take the letters in the 4 5 6 positions of the name, put that in a well column
    
  file_directory$field = substr(file_directory$file,10,11) %>% as.numeric()
  file_directory$field = file_directory$field +1 #nikon labels starting from 0. This makes it go from 1
  file_directory$channeltemporary = substr(file_directory$file,15,15) %>% as.numeric() 
  file_directory$channeltemporary = file_directory$channeltemporary +1 #so field starts from 1
  file_directory$channel = paste("CH_",file_directory$channeltemporary, sep = "")
  
  
  file_directory$file = NULL #don't need this column anymore, all key info split into other columns
  file_directory$channeltemporary= NULL
  
  file_directory = as_tibble(file_directory)
  
  if(is.null(lookup))
  {
    return(file_directory)
  }
  
  file_directory = file_directory %>% pivot_wider(names_from = channel, values_from = path)
  
  # Attach the lockup table data
  file_directory = right_join(file_directory, lookup, by = "well", relationship = "many-to-many")
  
  
  # Move the colnames where they need to be
  file_directory = file_directory %>%
    mutate(ch_dapi = case_when(ch_dapi == 1 ~ CH_1,
                               ch_dapi == 2 ~ CH_2,
                               ch_dapi == 3 ~ CH_3,
                               ch_dapi == 4 ~ CH_4)) %>%
    mutate(ch_antibody = case_when(ch_antibody == 1 ~ CH_1,
                                   ch_antibody == 2 ~ CH_2,
                                   ch_antibody == 3 ~ CH_3,
                                   ch_antibody == 4 ~ CH_4))
  
    if ("ch_actin" %in% names(file_directory)) {
      file_directory = file_directory %>% 
        mutate(ch_actin = case_when(ch_actin == 1 ~ CH_1,
                                  ch_actin == 2 ~ CH_2,
                                  ch_actin == 3 ~ CH_3,
                                  ch_actin == 4 ~ CH_4))
    }
  
  file_directory <- file_directory %>%
    select(-any_of(c("CH_1", "CH_2", "CH_3", "CH_4")))
  
  
  return(file_directory)
  
}




# Jayden's channels -------------------------------------------------------
