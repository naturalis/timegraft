# With the crown fixed at 1: how far do chronos node ages depart from a clock, and how long does it take?
suppressMessages({library(ape); library(parallel)})
out <- commandArgs(TRUE)[1]
tr <- read.delim(file.path(out, "trees.tsv"), quote = "", stringsAsFactors = FALSE)
FLOOR <- 1e-5
prep <- function(nwk) { t <- read.tree(text = nwk); t$edge.length[is.na(t$edge.length) | t$edge.length < FLOOR] <- FLOOR
  if (!is.binary(t)) t <- multi2di(t, random = FALSE); t$edge.length[t$edge.length < FLOOR] <- FLOOR; t }
set.seed(42)
small <- tr[tr$n >= 4 & tr$n <= 60, ]; small <- small[sample(nrow(small), 200), ]
big <- tr[tr$genus %in% c(tr$genus[which.min(abs(tr$n - 100))], tr$genus[which.min(abs(tr$n - 200))], tr$genus[which.min(abs(tr$n - 300))]), ]
one <- function(row, models) {
  t <- prep(row$newick)
  do.call(rbind, lapply(models, function(p) {
    set.seed(1); t0 <- proc.time()[3]
    ch <- tryCatch(suppressWarnings(chronos(t, lambda = as.numeric(p[2]), model = p[1], quiet = TRUE,
                   calibration = makeChronosCalib(t, node = "root", age.min = 1, age.max = 1))), error = function(e) NULL)
    ok <- !is.null(ch)
    data.frame(genus = row$genus, n = row$n, model = p[1], lambda = as.numeric(p[2]), ok = ok, secs = proc.time()[3] - t0,
               len = if (ok) sum(ch$edge.length) else NA, phiic = if (ok) attr(ch, "PHIIC")$PHIIC else NA,
               rate_cv = if (ok) sd(attr(ch, "rates")) / mean(attr(ch, "rates")) else NA,
               newick = if (ok) write.tree(ch) else NA)
  }))
}
models <- list(c("clock", 1), c("correlated", 1), c("relaxed", 1), c("correlated", 0.1), c("correlated", 10))
res <- mclapply(seq_len(nrow(small)), function(i) one(small[i, ], models), mc.cores = 2)
write.table(do.call(rbind, res), file.path(out, "C_fixed.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
resb <- mclapply(seq_len(nrow(big)), function(i) one(big[i, ], list(c("correlated", 1))), mc.cores = 2)
write.table(do.call(rbind, resb)[, 1:9], file.path(out, "C_big.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
