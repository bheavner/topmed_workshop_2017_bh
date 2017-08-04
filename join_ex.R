library(tidyverse)

parsed_snp_file <- "parsed_snp.tsv"
parsed_indel_file <- "parsed_indel.tsv"
unit_defs_file <- "feature_bounds_20170804.tsv"

snps <- read_tsv(parsed_snp_file, comment = "#")

indels <- read_tsv(parsed_indel_file, comment = "#")

unit_defs <- read_tsv(unit_defs_file, comment = "#", skip = 1)

unit_defs <- select(unit_defs, c(gene_id, agg_start, agg_end))

#this works - but how add aggregation unit names?
select(snps, c(chr, pos, ref, alt), num_range(unit_defs$agg_start, unit_defs$agg_end))

#there's probably a nice vectorized way to do this, but for demonstration purposes, this will work:
# make an empty tibble
foo <- tibble(group_id="", chromosome="", position="", ref="", alt="") %>%
  filter(length(group_id)>1)

# loop over unit defs
for (rowIndex in 1:nrow(unit_defs)) {

  # select snps and insert to foo
  snpsToAdd <- select(snps, c(chr, pos, ref, alt)) %>%
    dplyr::filter(between(pos, unit_defs[rowIndex, ]$agg_start, unit_defs[rowIndex, ]$agg_end)) %>%
    distinct() %>%
    mutate(group_id = unit_defs[rowIndex, ]$gene_id)
  
  if (nrow(snpsToAdd) > 0) {
    foo <- add_row(
      foo,
      group_id = snpsToAdd$group_id,
      chromosome = snpsToAdd$chr,
      position = snpsToAdd$pos,
      ref = snpsToAdd$ref,
      alt = snpsToAdd$alt
    )
  }
  
  # select indels and insert to foo
  toAdd <- select(indels, c(chr, pos, ref, alt)) %>%
    dplyr::filter(between(pos, unit_defs[rowIndex, ]$agg_start, unit_defs[rowIndex, ]$agg_end)) %>%
    distinct() %>%
    mutate(group_id = unit_defs[rowIndex, ]$gene_id)
  
  print(paste0("row: ", rowIndex, " snps to add: ", nrow(snpsToAdd), " indels to add: ", nrow(toAdd)))
  
  if (nrow(toAdd) > 0) {
    foo <- add_row(
      foo,
      group_id = toAdd$group_id,
      chromosome = toAdd$chr,
      position = toAdd$pos,
      ref = toAdd$ref,
      alt = toAdd$alt
    )
  }
}

glimpse(foo)

aggregated_variants <- distinct(foo)

#save it for the analysis pipelin!
save(aggregated_variants, file = "chr22_gene_aggregates.RDA")


# look at number of variants per aggregation unit:
counts <- aggregated_variants %>% group_by(group_id) %>% summarize(n())

range(counts$`n()`)
