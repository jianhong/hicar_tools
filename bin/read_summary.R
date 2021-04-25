#!/usr/bin/env Rscript
fs <- dir(".", "*.reads_stats.csv")
df <- do.call(rbind, lapply(fs, read.csv))
write.csv(df, "reads_summary.csv", row.names=FALSE)
