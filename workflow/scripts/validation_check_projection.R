#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
old <- read_tsv("results_project/admixture_cohort-ADSP-WGS-r4_ref-gnomad-hgdp-1kg.hg38.tsv")
new <- read_tsv("results/admixture_cohort-ADSP-WGS-r4_ref-gnomad-hgdp-1kg.hg38.tsv")
old_samp <- old |> filter(partition == "sample")
new_samp <- select(old_samp, IID) |> left_join(new, by = "IID")

message("Correlation tests:")
paste0("k", 1:12) |>
  map_dfr(
    ~ cor.test(pull(old_samp, .x), pull(new_samp, .x)) |>
      broom::tidy() |>
      mutate(k = .x)) |>
  select(k, everything())

message("Equality:")
tibble(k = paste0("k", 1:12)) |>
  mutate(
    equal = map_lgl(k, ~ all(pull(old_samp, .x) == pull(new_samp, .x))),
    equal_round_2 = map_lgl(k, ~ all(round(pull(old_samp, .x), digits = 2)
                                     == round(pull(new_samp, .x), digits = 2))),
    equal_round_1 = map_lgl(k, ~ all(round(pull(old_samp, .x), digits = 1)
                                     == round(pull(new_samp, .x), digits = 1))))

