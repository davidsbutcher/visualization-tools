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
# Note that this is a "best hits" waffle, i.e. unique proteins are NOT found
# in multiple fractions

shortnameToUse <- "gf05d"

# Packages ------------------------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(grDevices)
library(readxl)
library(waffle)
library(hrbrthemes)
library(ggthemes)
library(glue)
library(here)

# Load data for all shortnames ---------------------------------------------

# Filepath should be relative or absolute path to TDsummary file.

if (is.null(exclusionlist) == FALSE) {

TDsummary <- read_xlsx(TDsummarypath,
                       sheet = "summary")

TDsummary <- 
  TDsummary %>%
  dplyr::filter(is.na(`Short name`) != TRUE) %>%
  dplyr::filter(!(`Short name` %in% exclusionlist))    # Optional exclusion list

shortnamelist <- TDsummary$`Short name` %>% as.list

## Create tibbles for input to waffle function for proteins and proteoforms

waffledata_protein_all <- tibble(
  Localization = factor(
    c(
      rep("Cytosolic", length(TDsummary$`Short name`)),
      rep("Membrane", length(TDsummary$`Short name`)),
      rep("Periplasmic", length(TDsummary$`Short name`)),
      rep("NOTA", length(TDsummary$`Short name`))
    )
  ),
  Count = c(TDsummary$`Cytosolic Proteins`,
            TDsummary$`Membrane Proteins`,
            TDsummary$`Periplasmic Proteins`,
            TDsummary$`NOTA Proteins`),
  shortname = factor(rep(TDsummary$`Short name`, 4))
)

waffledata_proteoform_all <- tibble(
  Localization = factor(
    c(
      rep("Cytosolic", length(TDsummary$`Short name`)),
      rep("Membrane", length(TDsummary$`Short name`)),
      rep("Periplasmic", length(TDsummary$`Short name`)),
      rep("NOTA", length(TDsummary$`Short name`))
    )
  ),
  Count = c(TDsummary$`Cytosolic Proteoforms`,
            TDsummary$`Membrane Proteoforms`,
            TDsummary$`Periplasmic Proteoforms`,
            TDsummary$`NOTA Proteoforms`),
  shortname = factor(rep(TDsummary$`Short name`, 4))
)

} else {
  
  message("Exclusion list is set to NULL, skipping waffles based on shortnames")
  
}


# Load data by datafile ---------------------------------------------------

## Load protein counts split by fraction in TDreport

if (is.null(shortnameToUse) == FALSE) {
  
  TDsummary_byfraction <-
    read_xlsx(TDsummarypath,
              sheet = "summary_byfraction") %>% 
    drop_na() %>% 
    dplyr::filter(`Short name` == shortnameToUse)
  
  fractionlist <- 
    TDsummary_byfraction %>%
    .$fraction %>%
    as.list

## Create tibbles for input to waffle function for proteins
## split by fraction

waffledata_protein_all_byfraction <- 
  tibble(
    Localization = factor(
      c(
        rep("Cytosolic", length(fractionlist)),
        rep("Membrane", length(fractionlist)),
        rep("Periplasmic", length(fractionlist)),
        rep("NOTA", length(fractionlist))
      )
    ),
    Count = c(TDsummary_byfraction$`Cytosolic Proteins`,
              TDsummary_byfraction$`Membrane Proteins`,
              TDsummary_byfraction$`Periplasmic Proteins`,
              TDsummary_byfraction$`NOTA Proteins`),
    fraction = rep(TDsummary_byfraction$fraction, 4)
  ) %>%
  group_by(fraction, Localization) %>% 
  summarize(Count)

} else {
  
  message("Shortname to use is set to NULL, skipping waffles split by fraction")
  
}

# Make waffles for all shortnames -----------------------------------------

if (is.null(exclusionlist) == FALSE) {

## Protein Multi-waffle
  
output_waffle_multiple_protein <-
  waffledata_protein_all %>%
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
  waffledata_proteoform_all %>%
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


# Make waffle split by datafile/fraction ----------------------------------

if (is.null(shortnameToUse) == FALSE) {

output_waffle_protein_byfraction <-
  waffledata_protein_all_byfraction %>%
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


# Output, all shortname waffles -------------------------------------------

setwd(here())

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

# Output, waffle split by fraction ----------------------------------------

if (is_null(shortnameToUse) == FALSE) {
  
  message(glue("Saving waffle for {shortnameToUse} fractions to waffles directory"))
  
  #   png(filename = glue("output/waffles/{systime}_waffle_protein_byfraction_{shortnameToUse}.png"),
  #       width = 8, height = 5, units = "in",
  #       pointsize = 12, bg = "white", res = 600)
  #   print(output_waffle_multiple_protein2)
  #   dev.off()
  
  pdf(
    file = glue("output/waffles/{systime}_waffle_protein_byfraction_{shortnameToUse}.pdf"),
    width = 8, height = 5, bg = "transparent"
    )
  print(output_waffle_protein_byfraction)
  dev.off()
  
}