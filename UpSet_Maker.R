# This script is used for making UpSet plots from data in 
# EcoliMG1655_TDdatasummary.xlsx

# INITIAL PARAMETERS --------------------------------------------------------------------------

# TDsummarypath should be relative or absolute path to TDsummary file.
# This shouldn't need to be changed unless the MagLab Z drive is not mapped.
# Input should be xlsx file with sheet names formatted like 
# "{shortname}_fractions_protein" or "{shortname}_fractions_proteoform"
# which contains columns with UniProt IDs or proteoform record numbers

TDsummarypath <- "Z:/ICR/David Butcher/Data Summaries/EcoliMG1655_TDdatasummary.xlsx"

# List of shortnames of TD reports to get UpSet plots for

shortname <- c("gf03")

# Set the following parameters appropriately for the UpSet plots you
# want to make

make_protein_UpSet <- TRUE
make_proteoform_UpSet <- TRUE

# Install and Load Required Packages ----------------------------------------------------------

if(!require ("tidyverse")) install.packages("tidyverse")
if(!require ("magrittr")) install.packages("magrittr")
if(!require ("UpSetR")) install.packages("UpSetR")
if(!require ("grDevices")) install.packages("grDevices")
if(!require ("readxl")) install.packages("readxl")
if(!require ("glue")) install.packages("glue")
if(!require ("here")) install.packages("here")
if(!require ("svglite")) install.packages("svglite")

library(tidyverse)
library(magrittr)
library(UpSetR)
library(grDevices)
library(readxl)
library(glue)
library(here)
library(svglite)

# Functions -----------------------------------------------------------------------------------

saveProteinPlots <- function(plotobj, shortname, out_width = 8, out_height = 5,
                      dpi = 600, output_dir = getwd()) {

  png(filename = glue("{output_dir}{shortname}_UpSet_protein.png"),
      width = out_width, height = out_height, units = "in",
      pointsize = 12, bg = "white", res = dpi)
  print(plotobj)
  dev.off()
  
  pdf(file = glue("{output_dir}{shortname}_UpSet_protein.pdf"),
      width = out_width, height = out_height, bg = "transparent",
      useDingbats = FALSE)
  print(plotobj)
  dev.off()
  
  svg(file = glue("{output_dir}{shortname}_UpSet_protein.svg"),
          width = 8, height = 5,
          pointsize = 12, bg = "transparent")
  print(plotobj)
  dev.off()
  
}

saveProteoformPlots <- function(plotobj, shortname, width = 8, height = 5,
                             res = 600, output_dir = getwd()) {
  
  png(filename = glue("{output_dir}{shortname}_UpSet_proteoform.png"),
      width = width, height = height, units = "in",
      pointsize = 12, bg = "white", res = res)
  print(plotobj)
  dev.off()
  
  pdf(file = glue("{output_dir}{shortname}_UpSet_proteoform.pdf"),
      width = width, height = height, bg = "transparent",
      useDingbats = FALSE)
  print(plotobj)
  dev.off()
  
  svg(file = glue("{output_dir}{shortname}_UpSet_proteoform.svg"),
          width = 8, height = 5,
          pointsize = 12, bg = "transparent")
  print(plotobj)
  dev.off()
  
}

killNA <- function(list_in) {
  
  list_out <- list()
  
  for (i in seq_along(list_in)) {
  
   list_out[[i]] <- list_in[[i]][!is.na(list_in[[i]])]
  
  }
    
  names(list_out) <- names(list_in)
  
  return(list_out)
}
  

# Load Excel Sheets ---------------------------------------------------------------------------

# Proteins are listed by UniProt accession numbers,
# proteoforms are listed by proteoform record number.
# Sheets in TDsummary must be named correctly for this
# to work!

if (make_protein_UpSet == TRUE) {

protein_id <- shortname %>%
  as.list %>% 
  map(function(x) glue("{x}_fractions_protein")) %>%
  map(function(x) read_xlsx(TDsummarypath, sheet = x)) %>% 
  map(as.list) %>% 
  map(killNA)

names(protein_id) <- shortname

}

if (make_proteoform_UpSet == TRUE) {

proteoform_id <- shortname %>%
  as.list %>% 
  map(function(x) glue("{x}_fractions_proteoform")) %>%
  map(function(x) read_xlsx(TDsummarypath, sheet = x)) %>% 
  map(as.list) %>% 
  map(killNA)

names(proteoform_id) <- shortname

}

# UpSet plot, proteins by fraction ------------------------------------------------------------

if (make_protein_UpSet == TRUE) {

UpSet_protein <- 
  protein_id %>% 
  map(
    
    function(proteinlist)
    {
      
      upset(fromList(proteinlist),
            sets = rev(glue("Frac_0{1:length(proteinlist)}")),
            nintersects = NA,
            sets.x.label = "Total Protein IDs",
            keep.order = T,
            mainbar.y.label = "Unique Protein IDs in Intersection",
            text.scale =  c(1.5, 1.2, 1.5, 1.5, 1.2, 0.9),
            point.size = 2,
            line.size = 0.75,
            group.by = "degree")
      
    }
    
  )

}

# UpSet plot, proteoforms by fraction ---------------------------------------------------------

if (make_proteoform_UpSet == TRUE) {

UpSet_proteoform <- 
  proteoform_id %>% 
  map(
    
    function(proteoformlist)
    {
      
      upset(fromList(proteoformlist),
            sets = rev(glue("Frac_0{1:length(proteoformlist)}")),
            nintersects = NA,
            sets.x.label = "Total Proteoform IDs",
            keep.order = T,
            mainbar.y.label = "Unique Proteoform IDs in Intersection",
            text.scale =  c(1.5, 1.2, 1.5, 1.5, 1.2, 0.9),
            point.size = 2,
            line.size = 0.75,
            group.by = "degree")
      
    }
    
  )

}

# Output --------------------------------------------------------------------------------------

if (dir.exists("output/")== FALSE) {dir.create("output/")}
if (dir.exists("output/UpSet/")== FALSE) {dir.create("output/UpSet/")}

# To change output params of UpSet plots, check args for
# saveProteinPlots or saveProteoformPlots functions

shortnamelist <- shortname %>% as.list

if (make_protein_UpSet == TRUE) {
  
  message(glue("Saving protein UpSet plots to output/UpSet/")) 
  map2(UpSet_protein, shortnamelist, saveProteinPlots, output_dir = "output/UpSet/")
  
} else {
  
  message("Skipping protein UpSet plots!!")
  
}

if (make_proteoform_UpSet == TRUE) {
message(glue("Saving proteoform UpSet plots to output/UpSet/")) 
map2(UpSet_proteoform, shortnamelist, saveProteoformPlots, output_dir = "output/UpSet/")
} else {
  
  message("Skipping proteoform UpSet plots!!")
  
}