# Title: Driver enrichment analysis by age group in LUAD (AACR GENIE, v9.0)
# Description: Performs Fisher’s exact tests comparing driver mutation
# frequencies between patients <=45 and >45 years.
# Input: tab-delimited file with one column per driver and patient ages
# Output: odds ratios, confidence intervals, and p-values per driver

# Set useful options
options(scipen = 1000, stringsAsFactors = FALSE)

# ---- User inputs ----
input_file  <- "Landscape driver by age.tsv"
output_file <- "young_vs_old.le_45_years.driver_enrichments.sumstats.tsv"

# Read data as provided, filling NAs where necessary
x <- read.table(input_file, header = TRUE, fill = NA, sep = "\t")

# Melt dataframe such that one row = one patient
p <- as.data.frame(do.call("rbind", lapply(1:ncol(x), function(i){
  ages <- x[, i]
  ages <- ages[which(!is.na(ages))]
  data.frame("driver" = colnames(x)[i],
             "young" = ages <= 45)
})))

# Run one Fisher's exact test for each driver
res <- data.frame(do.call("rbind", lapply(colnames(x), function(driver){
  f.df <- data.frame("has_driver" = p$driver == driver,
                     "is_young" = p$young)
  unlist(fisher.test(table(f.df))[c("estimate", "conf.int", "p.value")])
})))

colnames(res) <- c("odds_ratio", "lower_95_ci", "upper_95_ci", "pvalue")
res$driver <- colnames(x)

# Write output
write.table(res[c("driver", setdiff(colnames(res), "driver"))], 
            output_file,
            sep = "\t", quote = FALSE, row.names = FALSE)
