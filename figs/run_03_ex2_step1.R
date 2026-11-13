# Step 1 of Exercise 2, notebook 03: CLR + Euclidean PCoA on genes
suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(microViz)
  library(phyloseq)
  library(patchwork)
})
source(here("R/count_tables.R"))

ps_genes <- as_counts(table = "rnaseq", min_count = 10, min_samples = 100)

# Rebuild the community column from notebook 02, Exercise 2: the Lactobacillus
# fraction of each sample's species reads, split at 0.5. (Saved to
# data/my_sample_groups.csv, which notebook 03's setup chunk reads.)
sp <- read_counts("instructional")
lacto_species <- sp$taxonomy$taxon[sp$taxonomy$Genus == "Lactobacillus"]
frac <- colSums(sp$species[rownames(sp$species) %in% lacto_species, , drop = FALSE]) /
  colSums(sp$species)
sample_groups <- tibble(
  sample_id      = names(frac),
  lacto_fraction = as.numeric(frac),
  community      = ifelse(lacto_fraction >= 0.5, "lactobacillus", "anaerobe")
)
write_csv(sample_groups, here("data/my_sample_groups.csv"))
rm(sp, frac, lacto_species); gc()

# attach community to the gene phyloseq object
sd <- as(sample_data(ps_genes), "data.frame") |>
  rownames_to_column("sid") |>
  left_join(select(sample_groups, sample_id, community), by = c("sid" = "sample_id")) |>
  column_to_rownames("sid")
sample_data(ps_genes) <- sample_data(sd)
rm(sd, sample_groups); gc()

cat("genes in ps_genes:", ntaxa(ps_genes), "| samples:", nsamples(ps_genes), "\n")
cat("ranks available:", paste(rank_names(ps_genes), collapse = ", "), "\n")

# ---- what pseudocount will CLR use? ----
# tax_agg sums genes sharing a symbol; CLR's zero fix (inside
# microbiome::transform, called by tax_transform) adds half of the smallest
# non-zero count of that aggregated table.
agg_mat <- as(otu_table(tax_agg(ps_genes, rank = "symbol")), "matrix")
cat("symbols after aggregation:", nrow(agg_mat), "\n")
cat("smallest non-zero aggregated count:", min(agg_mat[agg_mat > 0]), "\n")
cat("pseudocount CLR will add to zeros:", min(agg_mat[agg_mat > 0]) / 2, "\n")
cat("zeros in aggregated table:",
    sprintf("%d of %d entries (%.1f%%)", sum(agg_mat == 0), length(agg_mat),
            100 * sum(agg_mat == 0) / length(agg_mat)), "\n")
rm(agg_mat); gc()

# ---- the ordination ----
ord <- ps_genes |>
  tax_transform("clr", rank = "symbol") |>
  dist_calc("euclidean") |>
  ord_calc("PCoA")
gc()

eig <- ord@ord$eig
cat("axis 1 share:", round(100 * eig[1] / sum(eig), 1),
    "% | axis 2 share:", round(100 * eig[2] / sum(eig), 1), "%\n")

p1 <- ord_plot(ord, color = "community", size = 1.6, alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") + ggtitle("community")
p2 <- ord_plot(ord, color = "outcome", size = 1.6, alpha = 0.8) +
  scale_color_brewer(palette = "Set1") + ggtitle("outcome")
p3 <- ord_plot(ord, color = "rnaseq_batch", size = 1.6, alpha = 0.8) +
  scale_color_viridis_d() + ggtitle("rnaseq_batch")

p <- (p1 + p2 + p3) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(here("figs/03_ex2_step1_gene_pcoa.png"), p,
       width = 11, height = 4.4, dpi = 200, bg = "white")
cat("plot saved to figs/03_ex2_step1_gene_pcoa.png\n")

sd2 <- as(sample_data(ord), "data.frame")
for (v in c("community", "outcome", "rnaseq_batch")) {
  cat(sprintf("%-13s NA: %d\n", v, sum(is.na(sd2[[v]]))))
}
