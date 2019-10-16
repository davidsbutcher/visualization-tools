# This script will make a waffle plot for every shortname in
# EcoliMG1655_TDdatasummary.xlsx (aside from those in the exclusion list)
# or a waffle plot for a single shortname split by fractions (given by
# shortnameToUse)

# Initial Parameters --------------------------------------------------------------------------

# TDsummarypath should be relative or absolute path to TDsummary file.
# This shouldn't need to be changed unless the MagLab Z drive is not mapped.

TDsummarypath <- "Z:/ICR/David Butcher/Data Summaries/EcoliMG1655_TDdatasummary.xlsx"

# Shortnames to exclude when making a waffle including data from all shortnames
# in the TDsummary file. Set this to NULL to not make this waffle.

exclusionlist <- NULL

# Shortname to use when making waffle of a single sample set split by fraction
# Select only a single shortname! Set this to NULL to not make this waffle.

shortnameToUse <- "peppi04"

# Packages ------------------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(grDevices)
library(readxl)
library(waffle)
library(hrbrthemes)
library(ggthemes)
library(glue)

# Load Excel Data -----------------------------------------------------------------------------

# Filepath should be relative or absolute path to TDsummary file.

if (is.null(exclusionlist) == FALSE) {

TDsummary <- read_xlsx(TDsummarypath,
                       sheet = "summary")

TDsummary %<>%
  dplyr::filter(is.na(`Short name`) != TRUE) %>%
  dplyr::filter(!(`Short name` %in% exclusionlist))    # Optional exclusion list

shortnamelist <- TDsummary$`Short name` %>% as.list

## Create tibbles for input to waffle function for proteins and proteoforms

input_waffle_protein <- tibble(
  Localization = factor(
    c(
      rep("Cytosolic", length(TDsummary$`Short name`)),
      rep("Membrane", length(TDsummary$`Short name`)),
      rep("Both", length(TDsummary$`Short name`)),
      rep("None", length(TDsummary$`Short name`))
    )
  ),
  Count = c(TDsummary$`Cytosolic Proteins`,
            TDsummary$`Membrane Proteins`,
            TDsummary$`Cytosolic & Membrane Proteins`,
            TDsummary$`NoLoc Proteins`),
  shortname = factor(rep(TDsummary$`Short name`, 4))
)

input_waffle_proteoform <- tibble(
  Localization = factor(
    c(
      rep("Cytosolic", length(TDsummary$`Short name`)),
      rep("Membrane", length(TDsummary$`Short name`)),
      rep("Both", length(TDsummary$`Short name`)),
      rep("None", length(TDsummary$`Short name`))
    )
  ),
  Count = c(TDsummary$`Cytosolic Proteoforms`,
            TDsummary$`Membrane Proteoforms`,
            TDsummary$`Cytosolic & Membrane Proteoforms`,
            TDsummary$`NoLoc Proteoforms`),
  shortname = factor(rep(TDsummary$`Short name`, 4))
)

} else {
  
  message("Exclusion list is set to NULL, skipping waffles based on shortnames")
  
}

## Load protein counts split by datafile/fraction in TDreport

if (is.null(shortnameToUse) == FALSE) {

TDsummary_bydatafile <- read_xlsx(TDsummarypath,
                       sheet = "summary_bydatafile")

TDsummary_bydatafile %<>%
  dplyr::filter(is.na(`Short name`) != TRUE) %>%
  dplyr::filter(`Short name` == shortnameToUse)    # SELECT A SINGLE SHORTNAME!!

datafilelist <- TDsummary_bydatafile$raw_file_name %>% as.list

## Create tibbles for input to waffle function for proteins and proteoforms split by fraction

input_waffle_protein_byfraction <- tibble(
  Localization = factor(
    c(
      rep("Cytosolic", length(datafilelist)),
      rep("Membrane", length(datafilelist)),
      rep("Both", length(datafilelist)),
      rep("None", length(datafilelist))
    )
  ),
  Count = c(TDsummary_bydatafile$`Cytosolic Proteins`,
            TDsummary_bydatafile$`Membrane Proteins`,
            TDsummary_bydatafile$`Cytosolic & Membrane Proteins`,
            TDsummary_bydatafile$`NoLoc Proteins`),
  fraction = rep(TDsummary_bydatafile$fraction, 4)
)

} else {
  
  message("Shortname to use is set to NULL, skipping waffles split by fraction")
  
}

# Make Waffle Lists ---------------------------------------------------------------------------

input_protein_list <- list()

input_proteoform_list <- list()

input_protein_list2 <- list()


if (is.null(exclusionlist) == FALSE) {
  
  for (i in levels(input_waffle_protein$shortname)) {
    
    input_protein_list[[i]] <- input_waffle_protein %>% filter(shortname == i)
    
  }
  
  
  for (i in levels(input_waffle_proteoform$shortname)) {
    
    input_proteoform_list[[i]] <- input_waffle_proteoform %>% filter(shortname == i)
    
  }
  
}

if (is.null(shortnameToUse) == FALSE) {
  
  for (i in unique(input_waffle_protein_byfraction$fraction)) {
    
    input_protein_list2[[i]] <- input_waffle_protein_byfraction %>% filter(fraction == i)
    
  }
}

# Multiple Waffles ----------------------------------------------------------------------------

if (is.null(exclusionlist) == FALSE) {

output_waffle_multiple_protein <-
  input_waffle_protein %>%
  ggplot(aes(fill = Localization, values = Count)) +
  geom_waffle(color = "white", size = 0.35, n_rows = 10, flip = TRUE) +
  facet_wrap(~shortname, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(labels = function(x) x * 10, # this multiplier must be equal to n_rows above
                     expand = c(0,0)) +
  
  scale_fill_tableau(name=NULL) +
  coord_equal() +
  labs(
    x = "Short Name",
    y = "Protein Count"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line()) +
  guides(fill = guide_legend("Localization", reverse = TRUE))


## Proteoform Multi-waffle

output_waffle_multiple_proteoform <-
  input_waffle_proteoform %>%
  ggplot(aes(fill = Localization, values = Count)) +
  geom_waffle(color = "white", size = 0.35, n_rows = 10, flip = TRUE) +
  facet_wrap(~shortname, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(labels = function(x) x * 10, # this multiplier must be equal to n_rows above
                     expand = c(0,0)) +
  
  scale_fill_tableau(name=NULL) +
  coord_equal() +
  labs(
    x = "Short Name",
    y = "Proteoform Count"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line()) +
  guides(fill = guide_legend("Localization", reverse = TRUE))

}

## Waffle split by fraction

if (is.null(shortnameToUse) == FALSE) {

output_waffle_multiple_protein2 <-
  input_waffle_protein_byfraction %>%
  ggplot(aes(fill = Localization, values = Count)) +
  geom_waffle(color = "white", size = 0.35, n_rows = 10, flip = TRUE) +
  facet_wrap(~fraction, nrow = 1, strip.position = "bottom") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_continuous(labels = function(x) x * 10, # this multiplier must be equal to n_rows above
                     expand = c(0,0)) +
  
  scale_fill_tableau(name=NULL) +
  coord_equal() +
  labs(
    x = "Fraction Number",
    y = "Protein Count"
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(), axis.ticks.y = element_line()) +
  guides(fill = guide_legend("Localization", reverse = TRUE))

}

# Output --------------------------------------------------------------------------------------

# Check for the necessary output directories

if (dir.exists("output/")== FALSE) {dir.create("output/")}
if (dir.exists("output/waffles/")== FALSE) {dir.create("output/waffles/")}

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

if (is_null(exclusionlist) == FALSE) {
  
  message("Saving waffle for shortnames in TDsummary to output/waffles/ directory")
  
  #   png(filename = glue("output/waffles/{systime}_waffle_multiple_protein.png"),
  #       width = 8, height = 5, units = "in",
  #       pointsize = 12, bg = "white", res = 600)
  #   print(output_waffle_multiple_protein)
  #   dev.off()
  
  pdf(file = glue("output/waffles/{systime}_waffle_multiple_protein.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(output_waffle_multiple_protein)
  dev.off()
  
  
  # png(filename = glue("output/waffles/{systime}_waffle_multiple_proteoform.png"),
  #     width = 8, height = 5, units = "in",
  #     pointsize = 12, bg = "white", res = 600)
  # print(output_waffle_multiple_proteoform)
  # dev.off()
  
  pdf(file = glue("output/waffles/{systime}_waffle_multiple_proteoform.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(output_waffle_multiple_proteoform)
  dev.off()
  
}

# Waffles split by fraction -------------------------------------------------------------------

if (is_null(shortnameToUse) == FALSE) {
  
  message(glue("Saving waffle for {shortnameToUse} fractions to waffles directory"))
  
  #   png(filename = glue("output/waffles/{systime}_waffle_protein_byfraction_{shortnameToUse}.png"),
  #       width = 8, height = 5, units = "in",
  #       pointsize = 12, bg = "white", res = 600)
  #   print(output_waffle_multiple_protein2)
  #   dev.off()
  
  pdf(file = glue("output/waffles/{systime}_waffle_protein_byfraction_{shortnameToUse}.pdf"),
      width = 8, height = 5, bg = "transparent")
  print(output_waffle_multiple_protein2)
  dev.off()
  
}