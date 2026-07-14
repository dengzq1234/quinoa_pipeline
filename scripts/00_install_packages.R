install.packages("BiocManager", repos="https://cloud.r-project.org")
library(BiocManager)

BiocManager::install("biomaRt",        ask=FALSE, update=FALSE)
BiocManager::install("clusterProfiler", ask=FALSE, update=FALSE)
BiocManager::install("org.At.tair.db", ask=FALSE, update=FALSE)

install.packages("dplyr", repos="https://cloud.r-project.org")

message("=== Checking installations ===")
for (pkg in c("biomaRt", "clusterProfiler", "org.At.tair.db", "dplyr")) {
  status <- requireNamespace(pkg, quietly=TRUE)
  message(pkg, ": ", if (status) "OK" else "FAILED")
}
