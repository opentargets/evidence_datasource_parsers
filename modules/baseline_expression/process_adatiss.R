library(tidyverse)
library(arrow)

source('/home/alegbe/repos/evidence_datasource_parsers/modules/baseline_expression/AdaTiSS_fn.R')

dat.rna.path <- '/home/alegbe/results/unaggregated/dice/parquet/dice_baseline_expression'
dat.rna <- read_parquet(dat.rna.path)
X <- preproc.filter.fn(dat.rna, dat.type = "TPM or RPKM", proc.zero = "ceiled to 1", filter.col.prp = 1, exp.thres = 1)
p.dat <- read.csv(args[3])
tiss.abd <- tiss.abd.fn(X, p.dat)
result <- AdaTiSS(X, tiss.abd = tiss.abd, adjust = TRUE, adjust.opt = 0)
write.table(result$ada.s, file = args[4], sep = "\t", quote = FALSE, row.names = TRUE)
