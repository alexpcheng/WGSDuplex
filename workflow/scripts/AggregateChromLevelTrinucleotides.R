#!/usr/bin/env Rscript

library(data.table)

args = commandArgs(trailingOnly=TRUE)

output = args[1]
sample_id = args[2]
files = args[3:length(args)]
print(files)

df = rbindlist(lapply(X = files, FUN = fread))

df = aggregate(.~triN+variant, df[, c("triN", "variant", "count")], sum)
df$total_count = sum(df$count)
df$frequency = df$count / df$total_count

df$sample_id = sample_id

df = df[, c("sample_id", "triN", "variant", "frequency", "total_count", "count")]

fwrite(df, output, sep = '\t')
#print(df)

