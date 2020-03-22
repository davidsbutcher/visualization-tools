
library(tidyverse)
library(magrittr)
library(readxl)
library(glue)
library(here)
library(ggplot2)
library(furrr)

# Initial Params ----------------------------------------------------------

TDsummarypath <- "Z:/ICR/David Butcher/Data Summaries/EcoliMG1655_TDdatasummary.xlsx"

shortname <- c("peppi04b", "gf05c")

# Make workers for future_map

plan(multisession(workers = 2))

# Functions ---------------------------------------------------------------

killNA <- function(list_in) {
  
  list_out <- list()
  
  for (i in seq_along(list_in)) {
    
    list_out[[i]] <- list_in[[i]][!is.na(list_in[[i]])]
    
  }
  
  names(list_out) <- names(list_in)
  
  return(list_out)
}


# Analysis ----------------------------------------------------------------

TDsummary <- 
  read_xlsx(TDsummarypath, sheet = "summary") %>% 
  select(`Short name`, Lysate, `Fractionation Method`) %>% 
  mutate(
    Medium = case_when(str_detect(Lysate, "LB") ~ "LB",
                       str_detect(Lysate, "M9") ~ "M9",
                       TRUE ~ "NA")
  )

frac_methods <- 
  TDsummary %>% 
  pull("Fractionation Method") %>% 
  unique()

fraction_counts <- 
  shortname %>%
  as.list() %>% 
  map(function(x) glue("{x}_fractions_proteoform")) %>%
  future_map(function(x) read_xlsx(TDsummarypath, sheet = x)) %>% 
  map(as.list) %>% 
  map(killNA) %>% 
  map(unlist) %>% 
  map(unname) %>% 
  map(table) %>%
  map(enframe, name = "accession", value = "count") %>%
  set_names(shortname) %>% 
  map(group_by, count) %>% 
  map(summarize, proteoform_count = n()) %>% 
  map(rename, "fractions_found_in" = count) %>% 
  map2(names(.), ~mutate(.x, shortname = .y)) %>% 
  map(select, shortname, everything()) %>% 
  map(~mutate(.x, pform_frac = proteoform_count/sum(proteoform_count))) %>% 
  reduce(union_all) %>% 
  ungroup() %>% 
  left_join(TDsummary, by = c("shortname" = "Short name")) %>% 
  rename("frac_method" = `Fractionation Method`)



# Make Plots --------------------------------------------------------------

## Column plots - average of fractionation methods

column_plot_proteoform <- 
  fraction_counts %>%
  group_by(frac_method, Medium, fractions_found_in) %>% 
  summarize(
    ave_pform_frac = mean(pform_frac),
    stdev_pform_frac = sd(pform_frac)
  ) %>%
  ggplot(aes(fractions_found_in, ave_pform_frac, fill = frac_method)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(ymin = (ave_pform_frac - stdev_pform_frac/2),
        ymax = (ave_pform_frac + stdev_pform_frac/2)),
    width = 0.3, position = position_dodge(width = 0.9)
  ) +
  scale_x_continuous(
    breaks = seq(min(fraction_counts$fractions_found_in):max(fraction_counts$fractions_found_in))
  ) +
  labs(x = "# of Fractions Found In",
       y = "Fraction of Proteoform IDs",
       fill = "Fractionation\n Method") +
  theme(text = element_text(size = 24))

# Save Results ------------------------------------------------------------


