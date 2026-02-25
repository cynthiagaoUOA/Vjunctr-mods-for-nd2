#' Title
#'
#' @param mask
#'
#' @return
#' @export
#'
#' @examples
invert_mask = function(mask)
{
  (mask*-1)+1
}


#' Title
#'
#' @param image_data
#' @param mult
#'
#' @return
#' @export
#'
#' @examples
multiply = function(image_data, mult)
{
  image_data*mult
}


#' Title
#'
#' @param keyimage
#'
#' @return
#' @export
#'
#' @examples
tophat = function(keyimage, area = 50)
{

  ## Low-pass disc-shaped filter
  # f = makeBrush(100, shape='disc', step=FALSE)
  #
  # f = f/sum(f)
  y = gblur(keyimage, area)
  # plot(y*30)
  # plot(keyimage*30)

  k2 = (keyimage-y)
  # plot(k2*30)

  k2 = (k2>0)*k2

  k2 = as.Image(k2)

  return(k2)
}

# file_path = "E:\\Cynthia\\Stack test\\Cynthia[12001]\\water_plasmin_25_9[28601]\\2024-09-25T142749+1200[32701]\\2024-09-25T142749+1200[32701]"
# Operetta image import ---------------------------------------------------

import_operetta = function(file_path, lookup = NULL)
{

  # Grab the file list
  filelist = list.files(file_path, full.names = TRUE, pattern = ".tif")
  file_directory = data.frame(path = filelist)

  # Grab the metadata from the file and split it out
  file_directory$file = basename(file_path_sans_ext(file_directory$path))

  file_directory$well = paste(LETTERS[as.numeric(substr(file_directory$file,0,3))], substr(file_directory$file,5,6), sep = "")
  file_directory$field = substr(file_directory$file,8,8) %>% as.numeric()
  file_directory$channel = paste("CH_",substr(file_directory$file,16,18) %>% as.numeric(), sep = "")


  file_directory$plane = substr(file_directory$file,13,15) %>% as.numeric()


  file_directory$file = NULL

  file_directory = as.tibble(file_directory)

  if(is.null(lookup))
  {
    return(file_directory)
  }


  #QC
  file_directory
  unique(file_directory$well)

  file_directory = file_directory %>% pivot_wider(names_from = channel, values_from = path)

  # Attach the lockup table data
  file_directory = right_join(file_directory, lookup, by = "well")


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

  return(file_directory)

}



# imgdata = file_directory[1,]
# nuclear_disk = 10
# tophat_area = 10
# tophat_threshold = 0.005
# min_area = 20
# print_image = FALSE
# cropx = c(0,1)
# cropy = c(0,1)
# return_all = FALSE
# nuclear_area = 10
# nuclear_offset = 0.001



segment_and_quant = function(imgdata, nuclear_disk = 10, tophat_area = 10, tophat_threshold = 0.005, min_area = 20, nuclear_area = 40, nuclear_offset = 0.001, print_image = FALSE,
                             cropx = c(0,1), cropy = c(0,1), return_image = FALSE)
{

  # load images to RAM
  image = list()
  image[["stain"]] = readImage(imgdata$ch_antibody)
  image[["dapi"]] = readImage(imgdata$ch_dapi)
  image[["actin"]] = readImage(imgdata$ch_actin)


  # Crop a region if needs be
  xdirs = dim(image[["stain"]])[1] * cropx
  ydirs = dim(image[["stain"]])[2] * cropy


  image[["stain"]] = image[["stain"]][xdirs[1]:xdirs[2] , ydirs[1] : ydirs[2]]
  image[["dapi"]] = image[["dapi"]][xdirs[1]:xdirs[2] , ydirs[1] : ydirs[2]]
  image[["actin"]] = image[["actin"]][xdirs[1]:xdirs[2] , ydirs[1] : ydirs[2]]


  # View the target images
  # (image[["dapi"]]*10) %>% plot()
  # rgbImage(image[["stain"]]*10) %>% display()
  #
  # rgbImage(image[["stain"]]*10, image[["stain"]]*10, image[["stain"]]*10) %>% display()


  # Segment the nuclei

  kern = makeBrush(nuclear_disk, shape='disc')

  image[["mask_nucleus"]] = image[["dapi"]] %>% thresh(w = nuclear_area, h = nuclear_area, offset = nuclear_offset)

  # plot(image[["mask_nucleus"]] )

  image[["mask_nucleus_enlarged"]] = (image[["mask_nucleus"]]) %>% EBImage::dilate(kern)

  image[["mask_perinuclear"]] = (!(image[["mask_nucleus"]])) * image[["mask_nucleus_enlarged"]]


  image[["stain_no_nucleus"]] = (image[["mask_nucleus_enlarged"]] == 0) * image[["stain"]]

  # Make the cutout areas in question
  image[["stain_nuclear"]] = image[["stain"]] * image[["mask_nucleus"]]
  image[["stain_perinuclear"]] = (image[["stain"]] * image[["mask_nucleus_enlarged"]]) - (image[["stain"]] * image[["mask_nucleus"]])

  # paintObjects(image[["mask_nuclear_halo"]] %>% bwlabel(), toRGB(image[["dapi"]] *10), col = c("green"), opac = 0.8, thick = TRUE) %>% plot()

  # Segment the junctions

  image[["tophat"]] = image[["stain"]] %>% tophat(tophat_area)

  image[["tophat_no_nucleus"]] = image[["tophat"]] * (!image[["mask_nucleus_enlarged"]])

  # display(image[["tophat_no_nucleus"]]*10)

  # Threshold junctions
  image[["mask_stain_threshold"]] = image[["tophat_no_nucleus"]]>tophat_threshold


  image[["stain_threshold"]] = image[["stain_no_nucleus"]] * image[["mask_stain_threshold"]]
  image[["stain_background"]] = (image[["stain_no_nucleus"]] * !image[["mask_stain_threshold"]] * !image[["mask_nucleus_enlarged"]]) %>% EBImage::as.Image()


  # Remove any areas less than 100px

  image[["mask_stain_threshold_labeled"]] = image[["mask_stain_threshold"]] %>% bwlabel()
  area_sizes = imageData(image[["mask_stain_threshold_labeled"]]) %>% table() %>% as.data.frame() %>% mutate(ID = as.numeric(`.`)-1)


  areas_to_remove = subset(area_sizes, area_sizes$Freq<min_area)$ID

  image[["filtered_mask_stain_threshold_labeled"]]  = image[["mask_stain_threshold_labeled"]] %>% rmObjects(areas_to_remove)
  image[["mask_contiguous"]]  = image[["filtered_mask_stain_threshold_labeled"]]>0

  # Create a range of internal configurations useful data for later displays

  image[["mask_fleck"]] = image[["mask_stain_threshold"]] * !image[["mask_contiguous"]]

  image[["stain_contiguous"]] = (image[["stain"]] * image[["mask_contiguous"]])


  image[["stain_fleck"]] = (image[["stain"]] * image[["mask_stain_threshold"]]* !image[["mask_contiguous"]]) %>% as.Image()

  image[["stain_background"]] = (image[["stain"]]* !image[["mask_stain_threshold"]] %>% as.Image)  * (!image[["mask_nucleus"]] %>% as.Image())

  nd = imgdata

  nd$overall_actin= sum(image[["actin"]])
  nd$overall_dapi= sum(image[["dapi"]])

  nd$overall_stain= sum(image[["stain"]])
  nd$contiguous_area = sum(image[["mask_contiguous"]])
  nd$contiguous_fluorescence = sum(image[["stain_contiguous"]])
  nd$contiguous_fluorescence_per_area = nd$contiguous_fluorescence/nd$contiguous_area
  nd$nuclear_fluorescence = sum(image[["stain_nuclear"]])
  nd$perinuclear_fluorescence = sum(image[["stain_perinuclear"]])
  nd$fleck_fluorescence = sum(image[["stain_fleck"]])
  nd$background_fluorescence = sum(image[["stain_background"]])

  image[["stain_nuclear"]] %>% multiply(10) %>% display()



  if(return_image)
  {
    image[["quant"]] = nd
    return(image)
  }


  return(nd)


}






segment_and_quant_p = function(file_directory, nuclear_disk = 10, tophat_area = 10, tophat_threshold = 0.005, min_area = 20, nuclear_area = 40, nuclear_offset = 0.001, print_image = FALSE)
{


  library(doFuture)
  registerDoFuture()
  plan("multisession")

  library(progressr)
  handlers(global = TRUE)

  handlers(list(
    handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 100,
      complete = "+"
    )
  ))

  row_vector = c(1:nrow(file_directory))

  with_progress({

    p <- progressr::progressor(along = row_vector)

    output = foreach (kk =  row_vector) %dopar% {

      output_data = segment_and_quant(file_directory[kk,], nuclear_disk = nuclear_disk, tophat_area = tophat_area, tophat_threshold = tophat_threshold, min_area = min_area, nuclear_area = nuclear_area, print_image = FALSE)
      p(kk)
      return(output_data)
    }

    all_output = data.table::rbindlist(output) %>% as_tibble()

  })

  return(all_output)
}


# Image QC

segment_and_quant_i = function(
    file_directory,
    nuclear_disk = 10,
    tophat_area = 10,
    tophat_threshold = 0.005,
    min_area = 20,
    nuclear_area = 40,
    nuclear_offset = 0.001,
    print_image = FALSE){

  ui <- page_sidebar(
    title = "Vjunctur dashboard",
    sidebar = sidebar(
      title = "Histogram controls",
      selectInput(
        "output_type", "Select variable",
        c("rgb", "stain", "dapi", "actin", "segment","stain_perinuclear", "stain_nuclear", "stain_contiguous", "stain_background")
      ),
      sliderInput("imagechoice",
                  "Image to select:",
                  min = 1,
                  max = nrow(file_directory),
                  value = 1,
                  step = 1,
                  round = TRUE),
      numericInput("nuclear_disk",
                   "Size of nuclear disk",
                   value = nuclear_disk),
      numericInput("nuclear_area",
                   "nuclear_area",
                   value = nuclear_area),
      numericInput("nuclear_offset",
                   "nuclear_offset",
                   value = nuclear_offset),
      numericInput("tophat_area",
                   "tophat_area",
                   value = nuclear_disk),
      numericInput("tophat_threshold",
                   "tophat_threshold",
                   value = tophat_threshold),
      numericInput("min_area",
                   "min_area",
                   value = min_area),
      numericInput("brighten",
                   "brighten",
                   value = 1),
      sliderInput("cropy", "cropy",
                  max = 1, min = 0, value = c(0,1)),
      sliderInput("cropx", "cropx",
                  max = 1, min = 0, value = c(0,1))
    ),
    card(
      card_header("Images"),
      tableOutput("context"),
      displayOutput("widget")
    )
  )


  server <- function(input, output) {

    output$context = renderTable({
      calculate_image()[["quant"]] %>% select(-ch_dapi, -ch_actin, -ch_antibody)
    })



    calculate_image = reactive({
      image_stack = segment_and_quant(file_directory[input$imagechoice,],
                                      nuclear_disk = input$nuclear_disk,
                                      tophat_area = input$tophat_area,
                                      tophat_threshold = input$tophat_threshold,
                                      min_area = input$min_area,
                                      cropy = input$cropy,
                                      cropx = input$cropx,
                                      return_image = TRUE)

      print(paste("segment_and_quant(file_directory, nuclear_disk = ",input$nuclear_disk,", tophat_area = ",input$tophat_area,", tophat_threshold = ",input$tophat_threshold,", min_area = ",input$min_area,", nuclear_area = ",input$nuclear_area,", nuclear_offset = ",input$nuclear_offset,")"), sep = "")
      return(image_stack)
    })

    calculate_rgb = reactive({
      rgb = rgbImage(calculate_image()[["stain"]] %>% multiply(10*input$brighten),
                     calculate_image()[["actin"]] %>% multiply(10*input$brighten),
                     calculate_image()[["dapi"]] %>% multiply(10*input$brighten))
      return(display(rgb))
    })

    calculate_segment = reactive({
      rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(input$brighten))
      return(display(rgb))}
    )

    # calculate_segment = reactive({
    #   rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(1/0.2126*input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.2126*input$brighten) ,
    #                  calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(1/0.7152*input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.7152*input$brighten) ,
    #                  calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(0.5/0.0722*input$brighten))
    #   return(display(rgb))}
    # )


    # Generate the EBImage viewer here
    output$widget <- renderDisplay({

      if(input$output_type == "rgb")
      {
        return(calculate_rgb())
      }
      if(input$output_type == "segment")
      {
        return(calculate_segment())
      }


      display(calculate_image()[[input$output_type]] %>% multiply(10))
    })

  }

  # Run the application
  shinyApp(ui = ui, server = server)
}


# segment_and_quant_i(file_directory)














# I need to adjust this code so that it can handle not having actin
# and so that it can brighten the antibody stain


segment_and_quant_i = function(
    file_directory,
    nuclear_disk = 10,
    tophat_area = 10,
    tophat_threshold = 0.005,
    min_area = 20,
    nuclear_area = 40,
    nuclear_offset = 0.001,
    print_image = FALSE){
  
  ui <- page_sidebar(
    title = "Vjunctur dashboard",
    sidebar = sidebar(
      title = "Histogram controls",
      selectInput(
        "output_type", "Select variable",
        c("rgb", "stain", "dapi", "actin", "segment","stain_perinuclear", "stain_nuclear", "stain_contiguous", "stain_background")
      ),
      sliderInput("imagechoice",
                  "Image to select:",
                  min = 1,
                  max = nrow(file_directory),
                  value = 1,
                  step = 1,
                  round = TRUE),
      numericInput("nuclear_disk",
                   "Size of nuclear disk",
                   value = nuclear_disk),
      numericInput("nuclear_area",
                   "nuclear_area",
                   value = nuclear_area),
      numericInput("nuclear_offset",
                   "nuclear_offset",
                   value = nuclear_offset),
      numericInput("tophat_area",
                   "tophat_area",
                   value = nuclear_disk),
      numericInput("tophat_threshold",
                   "tophat_threshold",
                   value = tophat_threshold),
      numericInput("min_area",
                   "min_area",
                   value = min_area),
      numericInput("brighten",
                   "brighten",
                   value = 1),
      sliderInput("cropy", "cropy",
                  max = 1, min = 0, value = c(0,1)),
      sliderInput("cropx", "cropx",
                  max = 1, min = 0, value = c(0,1))
    ),
    card(
      card_header("Images"),
      tableOutput("context"),
      displayOutput("widget")
    )
  )
  
  
  server <- function(input, output) {
    
    output$context = renderTable({
      calculate_image()[["quant"]] %>% select(-ch_dapi, -ch_actin, -ch_antibody)
    })
    
    
    
    calculate_image = reactive({
      image_stack = segment_and_quant(file_directory[input$imagechoice,],
                                      nuclear_disk = input$nuclear_disk,
                                      tophat_area = input$tophat_area,
                                      tophat_threshold = input$tophat_threshold,
                                      min_area = input$min_area,
                                      cropy = input$cropy,
                                      cropx = input$cropx,
                                      return_image = TRUE)
      
      print(paste("segment_and_quant(file_directory, nuclear_disk = ",input$nuclear_disk,", tophat_area = ",input$tophat_area,", tophat_threshold = ",input$tophat_threshold,", min_area = ",input$min_area,", nuclear_area = ",input$nuclear_area,", nuclear_offset = ",input$nuclear_offset,")"), sep = "")
      return(image_stack)
    })
    
    calculate_rgb = reactive({
      rgb = rgbImage(calculate_image()[["stain"]] %>% multiply(10*input$brighten),
                     calculate_image()[["actin"]] %>% multiply(10*input$brighten),
                     calculate_image()[["dapi"]] %>% multiply(10*input$brighten))
      return(display(rgb))
    })
    
    calculate_segment = reactive({
      rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(10*input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(10*input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(10*input$brighten))
      return(display(rgb))}
    )
    
    # calculate_segment = reactive({
    #   rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(1/0.2126*input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.2126*input$brighten) ,
    #                  calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(1/0.7152*input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.7152*input$brighten) ,
    #                  calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(0.5/0.0722*input$brighten))
    #   return(display(rgb))}
    # )
    
    
    # Generate the EBImage viewer here
    output$widget <- renderDisplay({
      
      if(input$output_type == "rgb")
      {
        return(calculate_rgb())
      }
      if(input$output_type == "segment")
      {
        return(calculate_segment())
      }
      
      
      display(calculate_image()[[input$output_type]] %>% multiply(10*input$brighten))
    })
    
  }
  
  # Run the application
  shinyApp(ui = ui, server = server)
}



# seg and quant no actin --------------------------------------------------



segment_and_quant_i_noactin = function(
    file_directory,
    nuclear_disk = 10,
    tophat_area = 10,
    tophat_threshold = 0.005,
    min_area = 20,
    nuclear_area = 40,
    nuclear_offset = 0.001,
    print_image = FALSE){
  
  ui <- page_sidebar(
    title = "Vjunctur dashboard",
    sidebar = sidebar(
      title = "Histogram controls",
      selectInput(
        "output_type", "Select variable",
        c("rgb", "stain", "dapi", "segment","stain_perinuclear", "stain_nuclear", "stain_contiguous", "stain_background")
      ),
      sliderInput("imagechoice",
                  "Image to select:",
                  min = 1,
                  max = nrow(file_directory),
                  value = 1,
                  step = 1,
                  round = TRUE),
      numericInput("nuclear_disk",
                   "Size of nuclear disk",
                   value = nuclear_disk),
      numericInput("nuclear_area",
                   "nuclear_area",
                   value = nuclear_area),
      numericInput("nuclear_offset",
                   "nuclear_offset",
                   value = nuclear_offset),
      numericInput("tophat_area",
                   "tophat_area",
                   value = nuclear_disk),
      numericInput("tophat_threshold",
                   "tophat_threshold",
                   value = tophat_threshold),
      numericInput("min_area",
                   "min_area",
                   value = min_area),
      numericInput("brighten",
                   "brighten",
                   value = 1),
      sliderInput("cropy", "cropy",
                  max = 1, min = 0, value = c(0,1)),
      sliderInput("cropx", "cropx",
                  max = 1, min = 0, value = c(0,1))
    ),
    card(
      card_header("Images"),
      tableOutput("context"),
      displayOutput("widget")
    )
  )
  
  
  server <- function(input, output) {
    
    output$context = renderTable({
      calculate_image()[["quant"]] %>% select(-ch_dapi, -ch_antibody)
    })
    
    
    
    calculate_image = reactive({
      image_stack = segment_and_quant_noactin(file_directory[input$imagechoice,],
                                      nuclear_disk = input$nuclear_disk,
                                      tophat_area = input$tophat_area,
                                      tophat_threshold = input$tophat_threshold,
                                      min_area = input$min_area,
                                      cropy = input$cropy,
                                      cropx = input$cropx,
                                      return_image = TRUE)
      
      print(paste("segment_and_quant(file_directory, nuclear_disk = ",input$nuclear_disk,", tophat_area = ",input$tophat_area,", tophat_threshold = ",input$tophat_threshold,", min_area = ",input$min_area,", nuclear_area = ",input$nuclear_area,", nuclear_offset = ",input$nuclear_offset,")"), sep = "")
      return(image_stack)
    })
    
    calculate_rgb = reactive({
      rgb = rgbImage(calculate_image()[["stain"]] %>% multiply(10*input$brighten),

                     calculate_image()[["dapi"]] %>% multiply(10*input$brighten))
      return(display(rgb))
    })
    
    calculate_segment = reactive({
      rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(10*input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(10*input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5*input$brighten) ,
                     calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(10*input$brighten))
      return(display(rgb))}
    )
    
    # calculate_segment = reactive({
    #   rgb = rgbImage(calculate_image()[["stain_contiguous"]]%>% as.Image() %>% multiply(1/0.2126*input$brighten) + calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.2126*input$brighten) ,
    #                  calculate_image()[["stain_fleck"]] %>% as.Image() %>% multiply(1/0.7152*input$brighten)+ calculate_image()[["stain_nuclear"]]%>% as.Image() %>% multiply(0.5/0.7152*input$brighten) ,
    #                  calculate_image()[["stain_background"]]%>% as.Image() %>% multiply(0.5/0.0722*input$brighten))
    #   return(display(rgb))}
    # )
    
    
    # Generate the EBImage viewer here
    output$widget <- renderDisplay({
      
      if(input$output_type == "rgb")
      {
        return(calculate_rgb())
      }
      if(input$output_type == "segment")
      {
        return(calculate_segment())
      }
      
      
      display(calculate_image()[[input$output_type]] %>% multiply(10*input$brighten))
    })
    
  }
  
  # Run the application
  shinyApp(ui = ui, server = server)
}


# -----------------------------------------------


segment_and_quant_noactin = function(imgdata, nuclear_disk = 10, tophat_area = 10, tophat_threshold = 0.005, min_area = 20, nuclear_area = 40, nuclear_offset = 0.001, print_image = FALSE,
                             cropx = c(0,1), cropy = c(0,1), return_image = FALSE)
{
  
  # load images to RAM
  image = list()
  image[["stain"]] = readImage(imgdata$ch_antibody)
  image[["dapi"]] = readImage(imgdata$ch_dapi)

  
  
  # Crop a region if needs be
  xdirs = dim(image[["stain"]])[1] * cropx
  ydirs = dim(image[["stain"]])[2] * cropy
  
  
  image[["stain"]] = image[["stain"]][xdirs[1]:xdirs[2] , ydirs[1] : ydirs[2]]
  image[["dapi"]] = image[["dapi"]][xdirs[1]:xdirs[2] , ydirs[1] : ydirs[2]]

  
  
  # View the target images
  # (image[["dapi"]]*10) %>% plot()
  # rgbImage(image[["stain"]]*10) %>% display()
  #
  # rgbImage(image[["stain"]]*10, image[["stain"]]*10, image[["stain"]]*10) %>% display()
  
  
  # Segment the nuclei
  
  kern = makeBrush(nuclear_disk, shape='disc')
  
  image[["mask_nucleus"]] = image[["dapi"]] %>% thresh(w = nuclear_area, h = nuclear_area, offset = nuclear_offset)
  
  # plot(image[["mask_nucleus"]] )
  
  image[["mask_nucleus_enlarged"]] = (image[["mask_nucleus"]]) %>% EBImage::dilate(kern)
  
  image[["mask_perinuclear"]] = (!(image[["mask_nucleus"]])) * image[["mask_nucleus_enlarged"]]
  
  
  image[["stain_no_nucleus"]] = (image[["mask_nucleus_enlarged"]] == 0) * image[["stain"]]
  
  # Make the cutout areas in question
  image[["stain_nuclear"]] = image[["stain"]] * image[["mask_nucleus"]]
  image[["stain_perinuclear"]] = (image[["stain"]] * image[["mask_nucleus_enlarged"]]) - (image[["stain"]] * image[["mask_nucleus"]])
  
  # paintObjects(image[["mask_nuclear_halo"]] %>% bwlabel(), toRGB(image[["dapi"]] *10), col = c("green"), opac = 0.8, thick = TRUE) %>% plot()
  
  # Segment the junctions
  
  image[["tophat"]] = image[["stain"]] %>% tophat(tophat_area)
  
  image[["tophat_no_nucleus"]] = image[["tophat"]] * (!image[["mask_nucleus_enlarged"]])
  
  # display(image[["tophat_no_nucleus"]]*10)
  
  # Threshold junctions
  image[["mask_stain_threshold"]] = image[["tophat_no_nucleus"]]>tophat_threshold
  
  
  image[["stain_threshold"]] = image[["stain_no_nucleus"]] * image[["mask_stain_threshold"]]
  image[["stain_background"]] = (image[["stain_no_nucleus"]] * !image[["mask_stain_threshold"]] * !image[["mask_nucleus_enlarged"]]) %>% EBImage::as.Image()
  
  
  # Remove any areas less than 100px
  
  image[["mask_stain_threshold_labeled"]] = image[["mask_stain_threshold"]] %>% bwlabel()
  area_sizes = imageData(image[["mask_stain_threshold_labeled"]]) %>% table() %>% as.data.frame() %>% mutate(ID = as.numeric(`.`)-1)
  
  
  areas_to_remove = subset(area_sizes, area_sizes$Freq<min_area)$ID
  
  image[["filtered_mask_stain_threshold_labeled"]]  = image[["mask_stain_threshold_labeled"]] %>% rmObjects(areas_to_remove)
  image[["mask_contiguous"]]  = image[["filtered_mask_stain_threshold_labeled"]]>0
  
  # Create a range of internal configurations useful data for later displays
  
  image[["mask_fleck"]] = image[["mask_stain_threshold"]] * !image[["mask_contiguous"]]
  
  image[["stain_contiguous"]] = (image[["stain"]] * image[["mask_contiguous"]])
  
  
  image[["stain_fleck"]] = (image[["stain"]] * image[["mask_stain_threshold"]]* !image[["mask_contiguous"]]) %>% as.Image()
  
  image[["stain_background"]] = (image[["stain"]]* !image[["mask_stain_threshold"]] %>% as.Image)  * (!image[["mask_nucleus"]] %>% as.Image())
  
  nd = imgdata
  

  nd$overall_dapi= sum(image[["dapi"]])
  
  nd$overall_stain= sum(image[["stain"]])
  nd$contiguous_area = sum(image[["mask_contiguous"]])
  nd$contiguous_fluorescence = sum(image[["stain_contiguous"]])
  nd$contiguous_fluorescence_per_area = nd$contiguous_fluorescence/nd$contiguous_area
  nd$nuclear_fluorescence = sum(image[["stain_nuclear"]])
  nd$perinuclear_fluorescence = sum(image[["stain_perinuclear"]])
  nd$fleck_fluorescence = sum(image[["stain_fleck"]])
  nd$background_fluorescence = sum(image[["stain_background"]])
  
  image[["stain_nuclear"]] %>% multiply(10) %>% display()
  
  
  
  if(return_image)
  {
    image[["quant"]] = nd
    return(image)
  }
  
  
  return(nd)
  
  
}


segment_and_quant_p_noactin = function(file_directory, nuclear_disk = 10, tophat_area = 10, tophat_threshold = 0.005, min_area = 20, nuclear_area = 40, nuclear_offset = 0.001, print_image = FALSE)
{
  
  
  library(doFuture)
  registerDoFuture()
  plan("multisession")
  
  library(progressr)
  handlers(global = TRUE)
  
  handlers(list(
    handler_progress(
      format   = ":spin :current/:total (:message) [:bar] :percent in :elapsed ETA: :eta",
      width    = 100,
      complete = "+"
    )
  ))
  
  row_vector = c(1:nrow(file_directory))
  
  with_progress({
    
    p <- progressr::progressor(along = row_vector)
    
    output = foreach (kk =  row_vector) %dopar% {
      
      output_data = segment_and_quant_noactin(file_directory[kk,], nuclear_disk = nuclear_disk, tophat_area = tophat_area, tophat_threshold = tophat_threshold, min_area = min_area, nuclear_area = nuclear_area, print_image = FALSE)
      p(kk)
      return(output_data)
    }
    
    all_output = data.table::rbindlist(output) %>% as_tibble()
    
  })
  
  return(all_output)
}
