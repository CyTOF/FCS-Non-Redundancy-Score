### Opens all FCS files in directory as Flowset
### Arcsinh transforms according to data type (e.g. mass vs. flow)
### Calculates NRS for markers (non-redundancy score)
### Asks user which ones to plot

### Needs improving as to the extraction of marker names from FCS files?


#########################################################
### Installing and loading required packages
#########################################################

if (!require("svDialogs")) {
  install.packages("svDialogs", dependencies = TRUE)
  library(svDialogs)
}

if (!require("flowCore")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("flowCore")
    library(flowCore)
}


if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("reshape2")) {
  install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
}

if (!require("tcltk2")) {
  install.packages("tcltk2", dependencies = TRUE)
  library(tcltk2)
}



# Clear environment
rm(list = ls(all = TRUE))


# library(svDialogs)# Moved to top
# Get user input for file
testfile<-dlg_open()
# Convert to string value
testfile <- capture.output(testfile)[7]

if ((testfile)=="character(0)"){
  stop("File input cancelled")
}else{
  
  #Remove invalid characters from file input location
  testfile <- gsub("[\"]","",testfile)
  testfile<-substring (testfile,5)
  
  #Set file and directory
  filename <- basename (testfile)
  dir <- dirname (testfile)
  
  # Set working directory accoding to file chosen
  setwd(dir)
  
  # List FCS (and fcs!) files in this location
  filesToOpen <- unique(c(list.files(dir,pattern = ".FCS"),list.files(dir,pattern = ".fcs")))
  
  # Ask user if they want to open all files, or just some
  OpenAll <- askYesNo("Open all files in this directory?")
  
  if (OpenAll==FALSE){
    # Ask user which files to open
    filesToOpen <- tk_select.list(filesToOpen, multiple=TRUE,title="Select Files to open") 
  }
  
  # Load all into a flowset
  #library(flowCore) # Moved to top
  fcs_raw <- read.flowSet(filesToOpen, transformation = FALSE,
                          truncate_max_range = FALSE)
  
  # Create summary table
  panel_fcs <- pData(parameters(fcs_raw[[1]]))
  
  
  # Create list of filenames for the fcs raw data
  filenamelist <- rep(filesToOpen, fsApply(fcs_raw, nrow))
  
  # Determine whether data is CyTOF or Flow by presence of FSC
  # isflow will be 0 for a CyTOF or greater than 1 if flow
  isflow <- sum(grep("FSC",colnames(fcs_raw)))
  # Determine whether data is pre CyTOF 3 (Helios) by presence of "Cell_length", rather than "Event_length"
  isCyTOF2 <- sum(grep("Cell_length",colnames(fcs_raw)))
  # Determine if data is Helios
  isCyTOF3 <- sum(grep("Event_length",colnames(fcs_raw)))
  
  #asinh transform
  
  asinh_scale <- 5
  if (isflow>0){
    asinh_scale <- 150
  }
  
  fcs_raw <- fsApply(fcs_raw, function(x, cofactor = asinh_scale){
    #colnames(x) <- panel_fcs$desc
    expr <- exprs(x)
    expr <- asinh(expr / cofactor)
    exprs(x) <- expr
    x
  })
  
  
  ## Define a function that calculates the NRS (non-redundancy score) per sample
  NRS <- function(x, ncomp = 3){
    pr <- prcomp(x, center = TRUE, scale. = FALSE)
    score <- rowSums(outer(rep(1, ncol(x)),
                           pr$sdev[1:ncomp]^2) * abs(pr$rotation[,1:ncomp]))
    return(score)
  }
  
  ## Calculate the score
  nrs_sample <- fsApply(fcs_raw, NRS, use.exprs = TRUE)
  
  
  # Convert to data frame
  nrs_sample <- data.frame(nrs_sample)
  
  
  
  # Row names are the file names
  # rownames(nrs_sample) 
  
  # colnames are the marker names, including Time, Gaussian params etc.
  #colnames(nrs_sample)
  
  # If there are no desctiption parameters, use names instead
  #if (length(panel_fcs$desc[is.na(panel_fcs$desc)]) == length(panel_fcs$desc) == TRUE){
   # 
  #}
  
  # Otherwise, use descriptions
  if ((length(panel_fcs$desc[is.na(panel_fcs$desc)]) == length(panel_fcs$desc)) == FALSE){
    # Rename as per descriptions (instead of names)
    colnames(nrs_sample) <- panel_fcs$desc
    
    # Remove the blank columns (e.g. time, event length, gaussian)
    nrs_sample <- nrs_sample[,is.na(colnames(nrs_sample))==FALSE]
  }

  
  # Get list of mean nrs values per file
  nrs <- colMeans(nrs_sample, na.rm = TRUE)
  
  # Sort by decreasing NRS
  markers_ord <- names(sort(nrs, decreasing = TRUE))
  
  # Ask user which markers to plot
  markerlist <- tk_select.list(colnames(nrs_sample), multiple=TRUE,title="Select Markers to use (Cancel to use all).") 
  
  # If user cancels dialog box, use all markers.
  if(length(markerlist)==0 ){
    markerlist <- colnames(nrs_sample)
  }
  
  # Crop data to that
  nrs_sample <- nrs_sample[,markerlist]
  
  # Add Filename column
  nrs_sample$sample_id <- rownames(nrs_sample)
  
  # Melt into single table
  ggdf <- melt(nrs_sample, id.var = "sample_id",
               value.name = "nrs", variable.name = "antigen")
  
  # Markers as factors
  ggdf$antigen <- factor(ggdf$antigen, levels = markers_ord)
  

  ## Plot the NRS for ordered markers
  ggplot(ggdf, aes(x = antigen, y = nrs)) +
    geom_point(aes(color = sample_id), alpha = 0.9,
               position = position_jitter(width = 0.3, height = 0)) +
    geom_boxplot(outlier.color = NA, fill = NA) +
    stat_summary(fun.y = "mean", geom = "point", shape = 21, fill = "white") +
    theme_bw()+
    # comma separator for y axis
    #scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+
    # or log
    #scale_y_log10()+
    # Rotate X axis test
    theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
    theme(legend.position="top")
    

} # End of file cancel loop