# Does a range on the crown age give chronos an estimate? Rooted, undated genus trees of the Iberian run,
# crown bounded between half and twice the TimeTree crown age (and below the stem), three random starts.
# Usage: Rscript chronos_bounds.R OUT_DIR PROVENANCE MAX_TIPS   (MAX_TIPS was 300)
suppressMessages({library(ape); library(parallel)})
args <- commandArgs(TRUE)
out <- args[1]; prov <- read.delim(args[2], quote = "", stringsAsFactors = FALSE)
max_tips <- as.integer(args[3])
tr <- read.delim(file.path(out, "trees.tsv"), quote = "", stringsAsFactors = FALSE)
tr <- merge(tr, prov[, c("genus", "crown_age", "crown_source", "timetree_crown")], by = "genus")
FLOOR <- 1e-5  # about 1/150 of a substitution on 658 sites

prep <- function(nwk) {
  t <- read.tree(text = nwk)
  t$edge.length[is.na(t$edge.length) | t$edge.length < FLOOR] <- FLOOR
  if (!is.binary(t)) t <- multi2di(t, random = FALSE)
  t$edge.length[t$edge.length < FLOOR] <- FLOOR
  t
}
run <- function(t, model, lambda, amin, amax, seed) {
  set.seed(seed)
  cal <- makeChronosCalib(t, node = "root", age.min = amin, age.max = amax)
  t0 <- proc.time()[3]
  ch <- tryCatch(suppressWarnings(chronos(t, lambda = lambda, model = model, calibration = cal, quiet = TRUE)),
                 error = function(e) NULL)
  secs <- proc.time()[3] - t0
  if (is.null(ch)) return(list(ok = FALSE, secs = secs))
  bt <- branching.times(ch)
  list(ok = TRUE, secs = secs, root = max(bt), len = sum(ch$edge.length),
       ploglik = attr(ch, "ploglik"), phiic = attr(ch, "PHIIC")$PHIIC,
       rate_cv = { r <- attr(ch, "rates"); if (length(r) > 1) sd(r) / mean(r) else 0 })
}

# A. Does a range on the crown give an estimate? Bounds tt/2 .. 2*tt (below the stem), three starts.
A <- tr[tr$crown_source == "timetree" & tr$n >= 4 & tr$n <= max_tips, ]
resA <- mclapply(seq_len(nrow(A)), function(i) {
  t <- prep(A$newick[i]); tt <- A$timetree_crown[i]
  lo <- tt / 2; hi <- min(2 * tt, 0.9 * A$stem_age[i])
  do.call(rbind, lapply(c("correlated", "relaxed", "clock"), function(m)
    do.call(rbind, lapply(1:3, function(s) {
      r <- run(t, m, 1, lo, hi, s)
      data.frame(genus = A$genus[i], n = A$n[i], model = m, seed = s, lo = lo, hi = hi, tt = tt,
                 ok = r$ok, root = if (r$ok) r$root else NA, ploglik = if (r$ok) r$ploglik else NA)
    }))))
}, mc.cores = 2)
write.table(do.call(rbind, resA), file.path(out, "A_bounds.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
