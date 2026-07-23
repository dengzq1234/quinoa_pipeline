# Run once to install all packages required by the full pipeline.

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# ── Step 1: DESeq2 ─────────────────────────────────────────────────────────────
BiocManager::install("DESeq2",    ask = FALSE, update = FALSE)
BiocManager::install("apeglm",    ask = FALSE, update = FALSE)

# ── Step 2: Filter (optional) ──────────────────────────────────────────────────
if (!requireNamespace("rentrez", quietly = TRUE))
  install.packages("rentrez", repos = "https://cloud.r-project.org")
BiocManager::install("xml2", ask = FALSE, update = FALSE)

# ── Step 3: Annotation & GO ────────────────────────────────────────────────────
BiocManager::install("biomaRt",         ask = FALSE, update = FALSE)
BiocManager::install("clusterProfiler", ask = FALSE, update = FALSE)
BiocManager::install("org.At.tair.db",  ask = FALSE, update = FALSE)

# ── Shared utilities ───────────────────────────────────────────────────────────
for (pkg in c("dplyr", "ggplot2", "pheatmap", "tibble", "ggrepel", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

message("=== Checking installations ===")
for (pkg in c("DESeq2", "biomaRt", "clusterProfiler", "org.At.tair.db",
              "dplyr", "ggplot2", "pheatmap", "ggrepel", "stringr")) {
  message(pkg, ": ", if (requireNamespace(pkg, quietly = TRUE)) "OK" else "FAILED")
}
