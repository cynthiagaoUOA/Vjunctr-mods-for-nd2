# VJUNCTR for drugscreen paper testing. Cocktail A at timepoint 4

library(tidyverse)
library(tools)
source("96well_Ji_nikon_importfunction.R")
source("vjunctur_functions_v1.R")

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
library("EBImage")


# channel 0 is dapi, 1 is ve-cad, b-cat, then pecam
# ve is 546, bcat is 750

# increase all by one. dapi 1, vecad 2, b-cat 3, pecam 4



# import ------------------------------------------------------------------
T4_AJ_key <- tribble(~well, ~ch_dapi, ~ch_antibody, ~name_antibody, ~sample,
                     "B02", 1, 2, "VE-cad", "high melatonin",
                     "C02", 1, 2, "VE-cad", "high riluzole",
                     "D02", 1, 2, "VE-cad", "high cilostazol",
                     "E02", 1, 2, "VE-cad", "high ibuprofen",
                     "F02", 1, 2, "VE-cad", "high icatibant",
                     "G02", 1, 2, "VE-cad", "high pravastatin",
                     
                     "B02", 1, 3, "b-catenin", "high melatonin",
                     "C02", 1, 3, "b-catenin", "high riluzole",
                     "D02", 1, 3, "b-catenin", "high cilostazol",
                     "E02", 1, 3, "b-catenin", "high ibuprofen",
                     "F02", 1, 3, "b-catenin", "high icatibant",
                     "G02", 1, 3, "b-catenin", "high pravastatin",
                     
                     "B02", 1, 4, "PECAM", "high melatonin",
                     "C02", 1, 4, "PECAM", "high riluzole",
                     "D02", 1, 4, "PECAM", "high cilostazol",
                     "E02", 1, 4, "PECAM", "high ibuprofen",
                     "F02", 1, 4, "PECAM", "high icatibant",
                     "G02", 1, 4, "PECAM", "high pravastatin",
                     
                     ## 
                     "B04", 1, 2, "VE-cad", "low melatonin",
                     "C04", 1, 2, "VE-cad", "low riluzole",
                     "D04", 1, 2, "VE-cad", "low cilostazol",
                     "E04", 1, 2, "VE-cad", "low ibuprofen",
                     "F04", 1, 2, "VE-cad", "low icatibant",
                     "G04", 1, 2, "VE-cad", "low pravastatin",
                     
                     "B04", 1, 3, "b-catenin", "low melatonin",
                     "C04", 1, 3, "b-catenin", "low riluzole",
                     "D04", 1, 3, "b-catenin", "low cilostazol",
                     "E04", 1, 3, "b-catenin", "low ibuprofen",
                     "F04", 1, 3, "b-catenin", "low icatibant",
                     "G04", 1, 3, "b-catenin", "low pravastatin",
                     
                     "B04", 1, 4, "PECAM", "low melatonin",
                     "C04", 1, 4, "PECAM", "low riluzole",
                     "D04", 1, 4, "PECAM", "low cilostazol",
                     "E04", 1, 4, "PECAM", "low ibuprofen",
                     "F04", 1, 4, "PECAM", "low icatibant",
                     "G04", 1, 4, "PECAM", "low pravastatin",
                     
                     ##
                     "B06", 1, 2, "VE-cad", "high VPA",
                     "C06", 1, 2, "VE-cad", "high Rapamycin",
                     "D06", 1, 2, "VE-cad", "high Doxycycline",
                     "E06", 1, 2, "VE-cad", "high Fingolimod",
                     "F06", 1, 2, "VE-cad", "high Dipyridamole",
                     "G06", 1, 2, "VE-cad", "high Ticagrelor",
                     
                     "B06", 1, 3, "b-catenin", "high VPA",
                     "C06", 1, 3, "b-catenin", "high Rapamycin",
                     "D06", 1, 3, "b-catenin", "high Doxycycline",
                     "E06", 1, 3, "b-catenin", "high Fingolimod",
                     "F06", 1, 3, "b-catenin", "high Dipyridamole",
                     "G06", 1, 3, "b-catenin", "high Ticagrelor",
                     
                     "B06", 1, 4, "PECAM", "high VPA",
                     "C06", 1, 4, "PECAM", "high Rapamycin",
                     "D06", 1, 4, "PECAM", "high Doxycycline",
                     "E06", 1, 4, "PECAM", "high Fingolimod",
                     "F06", 1, 4, "PECAM", "high Dipyridamole",
                     "G06", 1, 4, "PECAM", "high Ticagrelor",
                     
                     ###
                     "B08", 1, 2, "VE-cad", "low VPA",
                     "C08", 1, 2, "VE-cad", "low Rapamycin",
                     "D08", 1, 2, "VE-cad", "low Doxycycline",
                     "E08", 1, 2, "VE-cad", "low Fingolimod",
                     "F08", 1, 2, "VE-cad", "low Dipyridamole",
                     "G08", 1, 2, "VE-cad", "low Ticagrelor",
                     
                     "B08", 1, 3, "b-catenin", "low VPA",
                     "C08", 1, 3, "b-catenin", "low Rapamycin",
                     "D08", 1, 3, "b-catenin", "low Doxycycline",
                     "E08", 1, 3, "b-catenin", "low Fingolimod",
                     "F08", 1, 3, "b-catenin", "low Dipyridamole",
                     "G08", 1, 3, "b-catenin", "low Ticagrelor",
                     
                     "B08", 1, 4, "PECAM", "low VPA",
                     "C08", 1, 4, "PECAM", "low Rapamycin",
                     "D08", 1, 4, "PECAM", "low Doxycycline",
                     "E08", 1, 4, "PECAM", "low Fingolimod",
                     "F08", 1, 4, "PECAM", "low Dipyridamole",
                     "G08", 1, 4, "PECAM", "low Ticagrelor",
                     
                     ## rest
                     "B10", 1, 2, "VE-cad", "Vehicle",
                     "C10", 1, 2, "VE-cad", "high imatinib",
                     "F10", 1, 2, "VE-cad", "high sapropterin",
                     
                     "B10", 1, 3, "b-catenin", "Vehicle",
                     "C10", 1, 3, "b-catenin", "high imatinib",
                     "F10", 1, 3, "b-catenin", "high sapropterin",
                     
                     "B10", 1, 4, "PECAM", "Vehicle",
                     "C10", 1, 4, "PECAM", "high imatinib",
                     "F10", 1, 4, "PECAM", "high sapropterin",
                  
                     ### 
                     "D10", 1, 2, "VE-cad", "low imatinib",
                     "G10", 1, 2, "VE-cad", "low sapropterin",
                     
                     "D10", 1, 3, "b-catenin", "low imatinib",
                     "G10", 1, 3, "b-catenin", "low sapropterin",
                     
                     "D10", 1, 4, "PECAM", "low imatinib",
                     "G10", 1, 4, "PECAM", "low sapropterin",
)
      

#import
T4AJs = CGimport_maxIPtif("drugscreenT4cocktailAvjunctr_singlechannelTIF", T4_AJ_key)
T4AJs$timepoint = "T4"

T24AJ = CGimport_maxIPtif("drugscreenT24cocktailAvjunctr_singlechannelTIF", T4_AJ_key)
T24AJ$timepoint = "T24"

# in cocktail B, 0 is claudin (1), 1 is zono occludin (2), 2 is dapi (3)
TJ_key <- tribble(~well, ~ch_dapi, ~ch_antibody, ~name_antibody, ~sample,
                     "B03", 3, 1, "Claudin-5", "high melatonin",
                     "C03", 3, 1, "Claudin-5", "high riluzole",
                     "D03", 3, 1, "Claudin-5", "high cilostazol",
                     "E03", 3, 1, "Claudin-5", "high ibuprofen",
                     "F03", 3, 1, "Claudin-5", "high icatibant",
                     "G03", 3, 1, "Claudin-5", "high pravastatin",
                     
                     "B03", 3, 2, "Zono occludin", "high melatonin",
                     "C03", 3, 2, "Zono occludin", "high riluzole",
                     "D03", 3, 2, "Zono occludin", "high cilostazol",
                     "E03", 3, 2, "Zono occludin", "high ibuprofen",
                     "F03", 3, 2, "Zono occludin", "high icatibant",
                     "G03", 3, 2, "Zono occludin", "high pravastatin",
                     
                     ## 
                     "B05", 3, 1, "Claudin-5", "low melatonin",
                     "C05", 3, 1, "Claudin-5", "low riluzole",
                     "D05", 3, 1, "Claudin-5", "low cilostazol",
                     "E05", 3, 1, "Claudin-5", "low ibuprofen",
                     "F05", 3, 1, "Claudin-5", "low icatibant",
                     "G05", 3, 1, "Claudin-5", "low pravastatin",
                     
                     "B05", 3, 2, "Zono occludin", "low melatonin",
                     "C05", 3, 2, "Zono occludin", "low riluzole",
                     "D05", 3, 2, "Zono occludin", "low cilostazol",
                     "E05", 3, 2, "Zono occludin", "low ibuprofen",
                     "F05", 3, 2, "Zono occludin", "low icatibant",
                     "G05", 3, 2, "Zono occludin", "low pravastatin",
                     
                     ##
                     "B07", 3, 1, "Claudin-5", "high VPA",
                     "C07", 3, 1, "Claudin-5", "high Rapamycin",
                     "D07", 3, 1, "Claudin-5", "high Doxycycline",
                     "E07", 3, 1, "Claudin-5", "high Fingolimod",
                     "F07", 3, 1, "Claudin-5", "high Dipyridamole",
                     "G07", 3, 1, "Claudin-5", "high Ticagrelor",
                    
                     "B07", 3, 2, "Zono occludin", "high VPA",
                     "C07", 3, 2, "Zono occludin", "high Rapamycin",
                     "D07", 3, 2, "Zono occludin", "high Doxycycline",
                     "E07", 3, 2, "Zono occludin", "high Fingolimod",
                     "F07", 3, 2, "Zono occludin", "high Dipyridamole",
                     "G07", 3, 2, "Zono occludin", "high Ticagrelor",
                     
                     
                     ###
                     "B09", 3, 1, "Claudin-5", "low VPA",
                     "C09", 3, 1, "Claudin-5", "low Rapamycin",
                     "D09", 3, 1, "Claudin-5", "low Doxycycline",
                     "E09", 3, 1, "Claudin-5", "low Fingolimod",
                     "F09", 3, 1, "Claudin-5", "low Dipyridamole",
                     "G09", 3, 1, "Claudin-5", "low Ticagrelor",
                     
                     "B09", 3, 2, "Zono occludin", "low VPA",
                     "C09", 3, 2, "Zono occludin", "low Rapamycin",
                     "D09", 3, 2, "Zono occludin", "low Doxycycline",
                     "E09", 3, 2, "Zono occludin", "low Fingolimod",
                     "F09", 3, 2, "Zono occludin", "low Dipyridamole",
                     "G09", 3, 2, "Zono occludin", "low Ticagrelor",
                     
                     
                     ## rest
                     "B11", 3, 1, "Claudin-5", "Vehicle",
                     "C11", 3, 1, "Claudin-5", "high imatinib",
                     "F11", 3, 1, "Claudin-5", "high sapropterin",
                     
                     "B11", 3, 2, "Zono occludin", "Vehicle",
                     "C11", 3, 2, "Zono occludin", "high imatinib",
                     "F11", 3, 2, "Zono occludin", "high sapropterin",
                     
                     ### 
                     "D11", 3, 1, "Claudin-5", "low imatinib",
                     "G11", 3, 1, "Claudin-5", "low sapropterin",
                    
                     "D11", 3, 2, "Zono occludin", "low imatinib",
                     "G11", 3, 2, "Zono occludin", "low sapropterin")
                     


#Tight junctions
T4TJs = CGimport_maxIPtif_threechannels("drugscreenT4cocktailBjunctr_singlechannelTIF", TJ_key)
T4TJs$timepoint = "T4"

T24TJ = CGimport_maxIPtif_threechannels("drugscreenT24cocktailBjunctr_singlechannelTIF", TJ_key)
T24TJ$timepoint = "T24"

TJs <- rbind(T4TJs, T24TJ)


# Making a plotting function once quant is performed ----------------------

normalise_summarise_plot <- function(dataset){
  # find means of vehicle
  water_means<- dataset  %>% filter(sample =="Vehicle") %>% 
    summarise(
      mean_cont_area = mean(contiguous_area),
      mean_cont_fluoro = mean(contiguous_fluorescence),
      mean_total= mean(overall_stain))
  
  # normalise all to vehicle by division
  normalised <- dataset %>% 
    mutate(norm_cont_area = contiguous_area/water_means$mean_cont_area,
           norm_cont_fluoro = contiguous_fluorescence/water_means$mean_cont_fluoro,
           norm_total= overall_stain/water_means$mean_total)
  
  # normalised as input for summarising
  
  summarised <- normalised %>% 
    group_by(sample) %>%  
    summarise(
      mean_cont_area = mean(norm_cont_area),
      sd_cont_area = sd(norm_cont_area),
      
      mean_total_fluoro = mean(norm_total),
      sd_total_fluoro = sd(norm_total),
      
      mean_cont_fluoro = mean(norm_cont_fluoro),
      sd_cont_fluoro = sd(norm_cont_fluoro))
  # no se bc only one run
  
  # cont area
  cont_area<- ggplot(summarised, aes(y=sample, x=mean_cont_area, group=sample))+
    geom_bar(stat="identity", fill= "lightgrey") + 
    geom_errorbar(
      aes(y= sample, xmax= mean_cont_area + sd_cont_area, xmin=mean_cont_area - sd_cont_area)) +
    geom_point(data=normalised,
               mapping = aes(
                 y=sample, 
                 x=norm_cont_area), 
               size = 1, colour= "blue", alpha=0.3) + theme_bw()+ 
    geom_vline(xintercept = 1, colour= "skyblue")+ labs(x="junctional area")
  
  # junctional fluoro
  cont_fluoro<- ggplot(summarised, aes(y=sample, x=mean_cont_fluoro, group=sample))+
    geom_bar(stat="identity", fill= "lightgrey") + 
    geom_errorbar(
      aes(y= sample, xmax= mean_cont_fluoro + sd_cont_fluoro, xmin=mean_cont_fluoro - sd_cont_fluoro)) +
    geom_point(data=normalised,
               mapping = aes(
                 y=sample, 
                 x=norm_cont_fluoro), 
               size = 1, colour= "blue", alpha=0.3) + theme_bw()+ 
    geom_vline(xintercept = 1, colour= "skyblue")+ labs(x="junctional fluorescence")
  
  # total fluoro
  total_fluoro<- ggplot(summarised, aes(y=sample, x=mean_total_fluoro, group=sample))+
    geom_bar(stat="identity", fill= "lightgrey") + 
    geom_errorbar(
      aes(y= sample, xmax= mean_total_fluoro + sd_total_fluoro, xmin=mean_total_fluoro - sd_total_fluoro)) +
    geom_point(data=normalised,
               mapping = aes(
                 y=sample, 
                 x=norm_total), 
               size = 1, colour= "blue",  alpha=0.3) + theme_bw()+ 
    geom_vline(xintercept = 1, colour= "skyblue")+ labs(x="total antibody fluorescence")
  
  library(patchwork)
  allplots<- cont_fluoro+ cont_area+ total_fluoro & labs(y=NULL)
  
  return(allplots)
}


# Plots AJs -------------------------------------------------------------------

# vecad T4
vecad <- T4AJs %>% filter(name_antibody== "VE-cad")

# segment_and_quant_i_noactin(vecad)

quant_vecad = segment_and_quant_p_noactin(
  vecad, nuclear_disk =  5 , tophat_area =  80 , tophat_threshold =  0.0001 ,min_area =  30 , nuclear_area =  2 , nuclear_offset =  0.001 
)

normalise_summarise_plot(quant_vecad)

# T24
vecad24 <- T24AJs %>% filter(name_antibody== "VE-cad")
# segment_and_quant_i_noactin(vecad24)
quant_vecad24 = segment_and_quant_p_noactin(
  vecad24, nuclear_disk =  5 , tophat_area =  80 , tophat_threshold =  0.0001 ,min_area =  30 , nuclear_area =  2 , nuclear_offset =  0.001 
)

normalise_summarise_plot(quant_vecad24)


## b-cat
bcat <- T4AJs %>% filter(name_antibody== "b-catenin")
segment_and_quant_i_noactin(bcat)

quant_bcat = segment_and_quant_p_noactin(
  bcat, nuclear_disk =  6 , tophat_area =  10 , tophat_threshold =  0.0002 ,min_area =  50 , nuclear_area =  40 , nuclear_offset =  0.001)


normalise_summarise_plot(quant_bcat)


bcat24 <- T24AJs %>% filter(name_antibody== "b-catenin")
segment_and_quant_i_noactin(bcat)

quant_bcat24 = segment_and_quant_p_noactin(
  bcat24, nuclear_disk =  6 , tophat_area =  10 , tophat_threshold =  0.0002 ,min_area =  50 , nuclear_area =  40 , nuclear_offset =  0.001)


normalise_summarise_plot(quant_bcat24)


## PECAM
pecam <- T4AJs %>% filter(name_antibody== "PECAM")
segment_and_quant_i_noactin(pecam)

quant_pecam = segment_and_quant_p_noactin(
  pecam, nuclear_disk =  10 , tophat_area =  8 , tophat_threshold =  0.0001 ,min_area =  50 , nuclear_area =  80 , nuclear_offset =  0.001)


normalise_summarise_plot(quant_pecam)

## T24
pecam24 <- T24AJs %>% filter(name_antibody== "PECAM")
segment_and_quant_i_noactin(pecam)

quant_pecam24 = segment_and_quant_p_noactin(
  pecam24, nuclear_disk =  10 , tophat_area =  8 , tophat_threshold =  0.0001 ,min_area =  50 , nuclear_area =  80 , nuclear_offset =  0.001)


normalise_summarise_plot(quant_pecam24)



# TJs ---------------------------------------------------------------------

claudin <- TJs %>% filter(name_antibody== "Claudin-5")

# segment_and_quant_i_noactin(claudin)

quant_claudin = segment_and_quant_p_noactin(
  claudin, nuclear_disk =  10 , tophat_area =  40, tophat_threshold =  0.002 ,min_area =  30, nuclear_area =  20, nuclear_offset =  0.001 
)

quant_claudin %>% filter(timepoint =="T24") %>% normalise_summarise_plot()

normalise_summarise_plot(quant_claudin)

## zono 
zono <- TJs %>% filter(name_antibody== "Zono occludin")
# segment_and_quant_i_noactin(zono)

quant_zono = segment_and_quant_p_noactin(
  zono, nuclear_disk =  15 , tophat_area =  50, tophat_threshold =  0.0002, min_area =  10, nuclear_area =  40, nuclear_offset =  0.001 
)

quant_claudin %>% filter(timepoint =="T24") %>% normalise_summarise_plot()


normalise_summarise_plot(quant_zono)


