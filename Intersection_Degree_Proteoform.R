
library(tidyverse)
library(magrittr)
library(readxl)
library(glue)
library(here)
library(ggplot2)
library(furrr)

# Initial Params ----------------------------------------------------------

TDsummarypath <- 
  "Z:/ICR/David Butcher/Data Summaries/EcoliMG1655/EcoliMG1655_TDdatasummary.xlsx"

shortname <- c("peppi04d")

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


# ggplot Themes -----------------------------------------------------------

column_plot_01 <- 
  list(
    scale_fill_viridis_d(
      option = "magma",
      begin = 0.25,
      end = 0.69),
    theme_minimal(),
    theme(
      panel.background = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      text = element_text(size = 16)
    )
  )

# Analysis ----------------------------------------------------------------

# Get fractionation method data from TDsummary and add column for
# for lysate medium

TDsummary <- 
  read_xlsx(TDsummarypath, sheet = "summary") %>% 
  select(`Short name`, Lysate, `Fractionation Method`) %>% 
  mutate(
    Medium = case_when(
      str_detect(Lysate, "LB") ~ "LB",
      str_detect(Lysate, "M9") ~ "M9",
      TRUE ~ "NA"
    )
  )

# Get list of fractionation methods from TDsummary

frac_methods <- 
  TDsummary %>% 
  pull("Fractionation Method") %>% 
  unique()

# 

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
  map(~mutate(.x, proteoform_frac = proteoform_count/sum(proteoform_count))) %>% 
  reduce(union_all) %>% 
  ungroup() %>% 
  left_join(TDsummary, by = c("shortname" = "Short name")) %>% 
  rename("frac_method" = `Fractionation Method`)

# Count number of GF and PEPPI proteins for use in column plot

gf_count <- 
  fraction_counts %>% 
  filter(frac_method == "GELFrEE") %>% 
  .$proteoform_count %>% 
  sum()

peppi_count <- 
  fraction_counts %>% 
  filter(frac_method == "PEPPI") %>% 
  .$proteoform_count %>% 
  sum()


# Make Plots --------------------------------------------------------------

## Summarize fraction counts for plotting

fraction_counts_summary <- 
  fraction_counts %>%
  group_by(frac_method, Medium, fractions_found_in) %>% 
  summarize(
    ave_proteoform_frac = mean(proteoform_frac),
    stdev_proteoform_frac = sd(proteoform_frac)
  ) 

## Column plots - average of fractionation methods

column_plot_proteoform <- 
  fraction_counts_summary %>% 
  ggplot(aes(fractions_found_in, ave_proteoform_frac*100, fill = frac_method)) +
  geom_col(position = "dodge") +
  geom_errorbar(
    aes(
      ymin = (ave_proteoform_frac*100 - (stdev_proteoform_frac/2)*100),
      ymax = (ave_proteoform_frac*100 + (stdev_proteoform_frac/2)*100)
    ),
    width = 0.3, position = position_dodge(width = 0.9)
  ) +
  scale_x_continuous(
    breaks = 
      seq(min(fraction_counts$fractions_found_in):max(fraction_counts$fractions_found_in))
  ) +
  scale_y_continuous(
    breaks = 
      seq(from = 0,
          to = plyr::round_any(max(fraction_counts$proteoform_frac*100), 10, ceiling),
          by = 5),
    expand = c(0,0)
  ) +
  labs(x = "Intersection Degree",
       y = "Percentage of Proteoform IDs",
       fill = "Fractionation\n Method") +
  # annotate("text", size = 4,
  #          x = max(fraction_counts$fractions_found_in - 0.5),
  #          y = max(fraction_counts$proteoform_frac*100*0.90),
  #          label = glue("GELFrEE: n = {gf_count}\n
  #                       PEPPI: n = {peppi_count}")
  # ) +
 column_plot_01

if (length(unique(fraction_counts_summary$frac_method)) == 1) {
  
  column_plot_proteoform <- 
    column_plot_proteoform +
    guides(fill = "none")
  
}


# Save Results ------------------------------------------------------------

if (dir.exists("output/")== FALSE) dir.create("output/")
if (dir.exists("output/int_degree/")== FALSE) dir.create("output/int_degree/")

systime <- format(Sys.time(), "%Y%m%d_%H%M%S")

ggsave(
  glue("output/int_degree/{systime}_intersection_degree_pform.png"),
  plot = column_plot_proteoform,
  width = 8,
  height = 5,
  dpi = 600
)

svg(
  file = glue("output/int_degree/{systime}_intersection_degree_pform.svg"), 
  width = 8,
  height = 5,
  pointsize = 12,
  bg = "transparent"
)
print(column_plot_proteoform)
dev.off()

pdf(
  file = glue("output/int_degree/{systime}_intersection_degree_pform.pdf"), 
  width = 8,
  height = 5,
  bg = "transparent",
  useDingbats = FALSE
)
print(column_plot_proteoform)
dev.off()
seq(min(fraction_counts$fractions_found_in):max(fraction_counts$fractions_found_in))

