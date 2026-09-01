if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("DESeq2",    ask = FALSE, update = FALSE)
BiocManager::install("apeglm",    ask = FALSE, update = FALSE)
BiocManager::install("DEGreport", ask = FALSE, update = FALSE)

for (pkg in c("dplyr", "ggplot2", "pheatmap", "tibble", "ggrepel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

cat("All packages installed.\n")
