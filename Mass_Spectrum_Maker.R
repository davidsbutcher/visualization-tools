# Packages ------------------------------------------------------------------------------------

library(here)
library(tictoc)
library(rawDiag)
library(magrittr)
library(tictoc)
library(glue)
library(fs)
library(tibble)
library(glue)
library(RSQLite)
library(purrr)
library(furrr)
library(readr)
library(openxlsx)
library(fs)
library(ggplot2)
library(stringr)
library(dplyr)

# Initialize Parameters ---------------------------------------------------

# Raw File Directory

rawFileDir <- 
  "Z:/ICR/David Butcher/21T data/"

#Raw File Name

rawFileName <- 
  "LCA_NT_20191113_PEPPI_EcoliRibosomes_F04_01.raw"

# Scan numbers to read

scanNum <- 
  c(453,454)

# m/z range to include. Set to NULL to use whole range

mzrange <- 
list(
  c(600, 2000)
)

# Output file size (width,height in inches)

outputSize <- c(10,6)

# Output file DPI

outputDPI <- 300

# Max number of future workers

maxWorkers <- 8

tic()

# Functions ---------------------------------------------------------------

kickout <- function(list) {
  
  # This function removes any element from the list of input files
  # (from root/input) which does not have one of the allowed
  # extensions or which has "deprecated"
  
  allowed_ext <- c("raw")
  
  for (i in rev(seq_along(list))) {
    
    if (!(tools::file_ext(list[[i]]) %in% allowed_ext)) {
      
      list[[i]] <- NULL 
      
    } else if (str_detect(list[[i]], fixed("deprecated", TRUE)) == TRUE) {
      
      list[[i]] <- NULL 
      
    }
  }
  
  return(list)
}


# Make future workers -----------------------------------------------------

if (length(scanNum) >= maxWorkers) {
  
  plan(multisession(workers = maxWorkers))
  
} else {
  
  plan(multisession(workers = length(scanNum)))
  
}

# Check filename and path -------------------------------------------------

rawFilesInDir <- 
  fs::dir_ls(
    rawFileDir, 
    recurse = TRUE,
    type = "file",
    regexp = c("[.]raw$")
  )

if (length(rawFilesInDir) == 0) {
  stop("No .raw files found in raw file directory")
}

if (rawFilesInDir %>% 
    str_detect(rawFileName) %>% 
    any() == FALSE) {
  stop("Input raw file not found in raw file directory")
}

rawFile <- 
  rawFilesInDir %>% 
  str_subset(rawFileName) %>% 
  as.list()


# Read raw file -----------------------------------------------------------

rawFileData <-
  map(rawFile, rawDiag::read.raw)

scanData <-
  map(
    scanNum,
    ~rawDiag::readScans(rawFile, .x)
  )

spectraList <-
  future_map(
    scanData,
    ~tibble(
      "mz" = .x[[1]]$mZ,
      "intensity" = .x[[1]]$intensity
    )
  ) %>%
  set_names(scanNum)

# ggplot Themes -----------------------------------------------------------

ms_theme <- 
  list(
    labs(
      x = "m/z",
      y = "Intensity"
    ),
    scale_y_continuous(
      expand = c(0,0),
      labels = scales::scientific_format(digits = 2)
    ),
    scale_x_continuous(breaks = scales::breaks_pretty(n = 6)),
    theme_minimal(),
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(size = 20)
    )
  )

ms_theme_02 <- 
  list(
    labs(
      x = "m/z",
      y = "Intensity"
    ),
    scale_x_continuous(breaks = scales::breaks_pretty(n = 6)),
    scale_y_continuous(expand = c(0,0)),
    theme_minimal(),
    theme(
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      text = element_text(size = 20)
    )
  )


# Make plots --------------------------------------------------------------

if (is.null(mzrange) == FALSE) {
  
  y_max <- 
    future_map2(
      spectraList,
      mzrange,
      ~.x %>% 
        filter(mz >= .y[[1]] & mz <= .y[[2]]) %>% 
        pull(intensity) %>% 
        max()
    )
  
  plotList <- 
    future_pmap(
      list(
        spectraList,
        mzrange,
        y_max
      ),
      ~ggplot(..1, aes(mz, intensity)) +
        geom_line(size = 0.5) +
        coord_cartesian(
          xlim = c(..2[[1]], ..2[[2]]),
          ylim = c(0, ..3 + ..3 * 0.05)
        ) +
        ms_theme_02
    )
  
} else {
  
  plotList <- 
    spectraList %>% 
    map(
      ~ggplot(.x, aes(mz, intensity)) +
        geom_line(size = 0.5) +
        ms_theme_02
      
    )
  
}

# Save plots --------------------------------------------------------------

plotFileNames <- 
  map2(
    rawFileName,
    scanNum,
    ~glue("output/mass_spectra/{fs::path_ext_remove(.x)}_scan{.y}.png")
  )

future_map2(
  plotFileNames,
  plotList,
  ~ggsave(
    .x,
    plot = .y,
    width = outputSize[1],
    height = outputSize[2],
    dpi = outputDPI,
    bg = "transparent"
  )
)

toc()