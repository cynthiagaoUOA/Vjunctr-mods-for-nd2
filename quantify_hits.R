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



# each antibody gets its own. Try making map like how I did with c5 and zono previously

T12_TJ_wellplate = tribble(~well, ~ch_dapi, ~ch_actin, ~ch_antibody, ~name_antibody, ~sample,
                        "B02", 1, 4, 2, "claudin-5", "low VPA vehicle",
                        "C02", 1, 4, 2, "claudin-5", "high VPA vehicle",
                        "D02", 1, 4, 2, "claudin-5", "low  VPA plasmin",
                        "E02", 1, 4, 2, "claudin-5", "high VPA plasmin",
                        "B02", 1, 4, 3, "zono occludin", "low VPA vehicle",
                        "C02", 1, 4, 3, "zono occludin", "high VPA vehicle",
                        "D02", 1, 4, 3, "zono occludin", "low  VPA plasmin",
                        "E02", 1, 4, 3, "zono occludin", "high VPA plasmin",
                        
                        "B04", 1, 4, 2, "claudin-5", "low LiCl vehicle",
                        "C04", 1, 4, 2, "claudin-5", "high LiCl vehicle",
                        "D04", 1, 4, 2, "claudin-5", "low  LiCl plasmin",
                        "E04", 1, 4, 2, "claudin-5", "high LiCl plasmin",
                        "B04", 1, 4, 3, "zono occludin", "low LiCl vehicle",
                        "C04", 1, 4, 3, "zono occludin", "high LiCl vehicle",
                        "D04", 1, 4, 3, "zono occludin", "low  LiCl plasmin",
                        "E04", 1, 4, 3, "zono occludin", "high LiCl plasmin",
                        
                        "B06", 1, 4, 2, "claudin-5", "low asta vehicle",
                        "C06", 1, 4, 2, "claudin-5", "high asta vehicle",
                        "D06", 1, 4, 2, "claudin-5", "low  asta plasmin",
                        "E06", 1, 4, 2, "claudin-5", "high asta plasmin",
                        "B06", 1, 4, 3, "zono occludin", "low asta vehicle",
                        "C06", 1, 4, 3, "zono occludin", "high asta vehicle",
                        "D06", 1, 4, 3, "zono occludin", "low  asta plasmin",
                        "E06", 1, 4, 3, "zono occludin", "high asta plasmin",
                        
                        "B10", 1,4,2, "claudin-5", "vehicle",
                        "B10", 1,4,3, "zono occludin", "vehicle",
                        
                        "D10", 1,4,2, "claudin-5", "plasmin",
                        "D10", 1,4,3, "zono occludin", "plasmin")
                        
# channel order verified on imageJ, goes: dapi cluadin zono actin
                        
T12hitsTJs = CGimport_maxIPtif("260213_hits_T12_TJs", 
                             T12_TJ_wellplate)
T12hitsTJs$timepoint = "T12"


# My test batch of two antibodies

# Quantify zono1 
zono <- T12hitsTJs %>% filter(name_antibody== "zono occludin")

segment_and_quant_i(zono)

quant_zono = segment_and_quant_p(
  zono, nuclear_disk =  10 , tophat_area =  100 , tophat_threshold =  0.0001 ,min_area =  10 , nuclear_area =  40 , nuclear_offset =  0.001 
)

ggplot(quant_zono, aes(x=sample, y=contiguous_fluorescence, group=sample))+
  geom_bar(stat="identity")


#claudin
c5 <- T12hitsTJs %>% filter(name_antibody== "claudin-5")

segment_and_quant_i(c5)

quant_claudin = segment_and_quant_p(
  c5, nuclear_disk =  5 , tophat_area =  15 , tophat_threshold =  0.001 ,min_area =  30 , nuclear_area =  5 , nuclear_offset =  0.001 
)


ggplot(quant_claudin, aes(x=sample, y=contiguous_fluorescence, group=sample))+
  geom_bar(stat="identity")



# --------- AJs

T12_AJ_wellplate = tribble(~well, ~ch_dapi, ~ch_antibody, ~name_antibody, ~sample,
                           "B03", 1,3, "ve-cadherin", "low VPA vehicle",
                           "C03", 1,3, "ve-cadherin", "high VPA vehicle",
                           "D03", 1,3, "ve-cadherin", "low  VPA plasmin",
                           "E03", 1,3, "ve-cadherin", "high VPA plasmin",
                           "B03", 1,2, "b-catenin", "low VPA vehicle",
                           "C03", 1,2, "b-catenin", "high VPA vehicle",
                           "D03", 1,2, "b-catenin", "low  VPA plasmin",
                           "E03", 1,2, "b-catenin", "high VPA plasmin",
                           "B03", 1,4, "pecam", "low VPA vehicle",
                           "C03", 1,4, "pecam", "high VPA vehicle",
                           "D03", 1,4, "pecam", "low  VPA plasmin",
                           "E03", 1,4, "pecam", "high VPA plasmin",
                           
                           "B05", 1,3, "ve-cadherin", "low LiCl vehicle",
                           "C05", 1,3, "ve-cadherin", "high LiCl vehicle",
                           "D05", 1,3, "ve-cadherin", "low  LiCl plasmin",
                           "E05", 1,3, "ve-cadherin", "high LiCl plasmin",
                           "B05", 1,2, "b-catenin", "low LiCl vehicle",
                           "C05", 1,2, "b-catenin", "high LiCl vehicle",
                           "D05", 1,2, "b-catenin", "low  LiCl plasmin",
                           "E05", 1,2, "b-catenin", "high LiCl plasmin",
                           "B03", 1,4, "pecam", "low LiCl vehicle",
                           "C03", 1,4, "pecam", "high LiCl vehicle",
                           "D03", 1,4, "pecam", "low  LiCl plasmin",
                           "E03", 1,4, "pecam", "high LiCl plasmin",
                           
                           "B07", 1,3, "ve-cadherin", "low asta vehicle",
                           "C07", 1,3, "ve-cadherin", "high asta vehicle",
                           "D07", 1,3, "ve-cadherin", "low  asta plasmin",
                           "E07", 1,3, "ve-cadherin", "high asta plasmin",
                           "B07", 1,2, "b-catenin", "low asta vehicle",
                           "C07", 1,2, "b-catenin", "high asta vehicle",
                           "D07", 1,2, "b-catenin", "low  asta plasmin",
                           "E07", 1,2, "b-catenin", "high asta plasmin",
                           "B03", 1,4, "pecam", "low asta vehicle",
                           "C03", 1,4, "pecam", "high asta vehicle",
                           "D03", 1,4, "pecam", "low  asta plasmin",
                           "E03", 1,4, "pecam", "high asta plasmin",
                           
                           "B11", 1,3, "ve-cadherin", "vehicle",
                           "B11", 1,2, "b-catenin", "vehicle",
                           "B11", 1,4, "pecam", "vehicle",
                           
                           "D11", 1,3, "ve-cadherin", "plasmin",
                           "D11", 1,2, "b-catenin", "plasmin",
                           "D11", 1,4, "pecam", "plasmin")

#segment and quant i breaks with AJs and not TJs. Because missing the actin. 



T12hitsAJs = CGimport_maxIPtif("260213_hits_T12_AJs", 
                             T12_AJ_wellplate)

T12hitsAJs$timepoint = "T12"

# VECAD
vecad <- T12hitsAJs %>% filter(name_antibody== "ve-cadherin")

segment_and_quant_i_noactin(vecad)



quant_vecad = segment_and_quant_p_noactin(
  vecad, nuclear_disk =  1 , tophat_area =  50 , tophat_threshold =  0.001 ,min_area =  20 , nuclear_area =  40 , nuclear_offset =  0.001 
)

ggplot(quant_vecad, aes(x=sample, y=contiguous_fluorescence, group=sample))+
  geom_bar(stat="identity")

vecad <- T12hitsAJs %>% filter(name_antibody== "ve-cadherin")

segment_and_quant_i_noactin(vecad)

# b-catenin

bcat <- T12hitsAJs %>% filter(name_antibody== "b-catenin")

segment_and_quant_i_noactin(bcat)

quant_bcat = segment_and_quant_p_noactin(
  bcat, nuclear_disk =  1 , tophat_area =  10 , tophat_threshold =  0.002 ,min_area =  20 , nuclear_area =  20 , nuclear_offset =  0.001 
)



ggplot(quant_bcat, aes(x=sample, y=contiguous_fluorescence, group=sample))+
  geom_bar(stat="identity")


