library(readr)
library(tidyr)
suppressPackageStartupMessages(library(dplyr))
library(purrr)

input_metadata <-
  "Downloads/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz"
output_table <- "reference_proc/hgdp_1kg.popdata.tsv.gz"
output_smap <- "reference_proc/hgdp_1kg.smap"

input_metadata <- snakemake@input[["metadata"]]
output_table <- snakemake@output[["table"]]
output_smap <- snakemake@output[["smap"]]

parse_gnomad <- Vectorize(function(x) {
    # if the input is NA, return NA
    if (is.na(x)) {
        return(NA)
    # otherwise, parse the JSON
    } else {
        return(jsonlite::fromJSON(x))
    }
}, SIMPLIFY = FALSE, USE.NAMES = FALSE)

poptab <- input_metadata |>
  read_tsv() |>
  select(s, gnomad_population_inference, hgdp_tgp_meta, high_quality) |>
  mutate(popinf = map_vec(gnomad_population_inference, parse_gnomad),
         metadata =  map_vec(hgdp_tgp_meta, parse_gnomad)) |>
  filter(!is.na(popinf) & !is.na(high_quality), high_quality) |>
  select(-high_quality, -gnomad_population_inference, -hgdp_tgp_meta) |>
  unnest_wider(all_of(c("popinf", "metadata"))) |>
  mutate(pop_proc = case_when(
           pop %in% c("fin", "nfe") ~ "EUR",
           genetic_region == "OCE" ~ "OCE",
           TRUE ~ toupper(pop)),
         spop_checked = case_when(
           genetic_region == "CSA" & pop_proc == "SAS" ~ "SAS",
           genetic_region == pop_proc ~ pop_proc,
           TRUE ~ NA_character_)) |>
  select(-hgdp_technical_meta, -global_pca_scores, -subcontinental_pca) |>
  rename(pc = pca_scores, pop_orig = pop, spop = pop_proc, pop = population,
         IID = s) |>
  unnest_wider(pc, names_sep = "") |>
  select(IID, pop, spop, everything()) |>
  write_tsv(output_table)

poptab |>
  select("#Sample" = "IID", "Panel" = "spop") |>
  write_tsv(output_smap)
