suppressPackageStartupMessages(library(dplyr))
library(readr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)

## Input and output files

if (exists("snakemake")) {
  in_q_ref <- snakemake@input[["q_ref"]]
  in_fam_ref <- snakemake@input[["fam_ref"]]
  in_q_samp <- snakemake@input[["q_samp"]]
  in_fam_samp <- snakemake@input[["fam_samp"]]
  in_pops <- snakemake@input[["pops"]]
  out_anc <- snakemake@output[["anc"]]
  out_rda <- snakemake@output[["rda"]]
} else {
  setwd("/sc/arion/projects/load/users/fultob01/projects/nadmixture")
  in_q <- "gnomad.genomes.v3.1.2.hgdp_tgp.pruned.12.Q"
  in_fam <- "filtered_gnomad_1kg_hgdp/gnomad.genomes.v3.1.2.hgdp_tgp.pruned.fam"
  in_pops <- "outputs/hgdp_1kg.popdata.tsv"
  out_smap <- "reference_proc/hgdp_1kg.filt.smap"
  out_rda <- "reference_proc/hgdp_1kg.admixture.rda"
  out_samplist <- "reference_proc/hgdp_1kg.filt.samplist"
}

# Fam and popfiles
## ======================================##
message("Reading pop file \n")
pops <- in_pops |>
  read_tsv(col_types = cols(.default = "c"))

message("Reading fam files \n")
read_fam <- function(in_fam) {
  in_fam |>
    read_table(col_names = c("IID"), col_types = "-c----") |>
    mutate(order = row_number())
}

famfile_ref <- read_fam(in_fam_ref) |>
  mutate(partition = "reference") |>
  left_join(pops, by = c("IID"))

famfile_samp <- read_fam(in_fam_samp) |>
  mutate(partition = "sample")

famfile <- bind_rows(famfile_ref, famfile_samp)

# Interpreting unsupervised admixture output #
## ======================================##
message("Reading unsupervised admixture output \n")

read_q <- function(in_q, fam) {
  in_q |>
    read_table(col_names = FALSE, col_types = cols(.default = "d")) |>
    bind_cols(fam) |>
    rename_with(~ str_replace(.x, "^X", "k"))
}

tbl_admix_ref <- read_q(in_q_ref, famfile_ref)
tbl_admix_samp <- read_q(in_q_samp, famfile_samp)
tbl_admix <- bind_rows(tbl_admix_ref, tbl_admix_samp)

# Determining cluster labels

cluster_cols <- names(tbl_admix)[str_detect(names(tbl_admix), "^k\\d+$")]

assign_labels <- function(tbl_admix) {
  if ("spop_checked" %in% colnames(tbl_admix)) {
    assign_admix_raw <- tbl_admix |>
      select(any_of(c("FID", "IID")),
        spop = spop_checked, matches("^k\\d+$")) |>
      filter(spop != "MID") |> # remove middle eastern from assignment
      group_by(spop) |>
      summarise(across(where(is.numeric), mean)) |>
      filter(!is.na(spop))
  } else {
    assign_admix_raw <- tbl_admix |>
      group_by(spop) |>
      summarise(across(where(is.numeric), mean)) |>
      filter(!is.na(spop))
  }

  assign_admix_mat <- assign_admix_raw |>
    column_to_rownames(var = "spop") |>
    as.matrix()

  assign_admix <- assign_admix_mat |>
    t() |>
    as.data.frame() |>
    (\(.) mutate(., anc = colnames(.)[apply(., 1, which.max)]))() |>
    as_tibble(rownames = "cluster") |>
    rowwise() |>
    mutate(maxval = max(c_across(where(is.numeric)))) |>
    group_by(anc) |>
    arrange(-maxval) |>
    mutate(n = n(),
           cname = ifelse(n > 1, paste(anc, row_number(), sep = "_"), anc)) |>
    ungroup() |>
    select(-maxval, -n) |>
    arrange(cname) |>
    select(cname, cluster, anc, everything())

  assign_cname_vec <- pull(assign_admix, cname, cluster)

  heatmap_names <- colnames(assign_admix_mat) |>
    (\(x) sprintf("%s (%s)", assign_cname_vec[x], x))()

  heatmap_mat <- assign_admix_mat
  colnames(heatmap_mat) <- heatmap_names

  return(list(assign = assign_admix,
              heatmap = heatmap_mat,
              assign_cname_vec = assign_cname_vec))
}

admix_labs <- assign_labels(tbl_admix)

assign_admix <- admix_labs$assign
heatmap_mat <- admix_labs$heatmap
assign_cname_vec <- admix_labs$assign_cname_vec
assign_super_vec <- pull(assign_admix, anc, cluster)
assign_cname_vec_i <- pull(assign_admix, cluster, cname)
superpops <- rownames(heatmap_mat)

rm(admix_labs, assign_labels)

# Assign individuals

collapse_superpop <- function(df, sp) {
  # Add overall proportion of each superpop, collapsing clusters
  get_clusters <- \(spop) names(assign_super_vec[assign_super_vec == spop])
  mutate(df, !!sp := rowSums(across(all_of(get_clusters(sp)))))
}

tbl_admix_collapsed <- tbl_admix
for (sp in superpops) {
  tbl_admix_collapsed <- collapse_superpop(tbl_admix_collapsed, sp)
}

tbl_admix_inf <- tbl_admix_collapsed |>
  rowwise(IID) |>
  mutate(maxval = max(c_across(all_of(cluster_cols))),
         matchval = which.max(c_across(all_of(cluster_cols))),
         max_spop_prop = max(c_across(all_of(superpops)))) |>
  ungroup() |>
  (\(.) mutate(.,
    maxclust = colnames(.)[max.col(select(., matches("^k\\d+$")))],
    "Maximum Cluster" = unname(assign_cname_vec[maxclust]),
    admixture_super_pop_max = map_chr(maxclust, \(x) assign_super_vec[[x]]),
    admixture_cluster_max = map_chr(maxclust, \(x) assign_cname_vec[[x]])))() |>
  arrange(spop, admixture_super_pop_max, matchval, -maxval) |>
  mutate(
    pop = forcats::fct_inorder(pop),
    spop = forcats::fct_inorder(spop),
    admixture_super_pop_max = factor(
      admixture_super_pop_max, levels = levels(spop))) |>
  arrange(matchval, -maxval) |>
  mutate(IID = forcats::fct_inorder(IID))

out_admix <- tbl_admix_inf |>
  select(-maxval, -matchval, -maxclust) |>
  select(any_of(c("FID", "IID")),
    all_of(superpops), matches("^k\\d+$"),
    everything()) |>
  filter(!is.na(pop) | partition == "sample")

# Filter samples

out_admix |>
  write_tsv(out_anc)

save(out_admix, heatmap_mat, assign_cname_vec, assign_cname_vec_i,
  superpops, cluster_cols, file = out_rda)
