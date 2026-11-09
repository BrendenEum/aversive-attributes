# Sample the two by-study attention-weight models (Narrow, Wide; the main
# model with every coefficient crossed with study; an Extra Analyses
# robustness check) outside of RStudio, so that a closed session can't kill a
# run that takes many hours. From the project root, at low priority so that
# foreground work keeps the CPU when it needs it (four chains, one core each;
# the two fits run one after the other):
#
#   nice -n 10 Rscript code/fitting/sample_bystudy_attention_weights.R > results/logs/sampling_utility_bystudy.log 2>&1 &
#
# With a range as the argument it samples only that fit:
#
#   nice -n 10 Rscript code/fitting/sample_bystudy_attention_weights.R Narrow > results/logs/sampling_bystudy_narrow.log 2>&1 &
#
# It evaluates the notebook's own chunks (libraries, data, design matrices,
# fit), so the caches it leaves in notebooks/brm_fits/ are exactly what
# Main_Analysis.Rmd loads afterwards (with run_extra_analyses <- TRUE). brm()
# skips any fit whose cache already exists, so the script can be re-run after
# an interruption.

rmd <- here::here("notebooks/Main_Analysis.Rmd")

# Split the notebook into its code chunks, keyed by chunk label.
lines  <- readLines(rmd)
starts <- grep("^```\\{r [^,}]+", lines)
ends   <- grep("^```$", lines)
chunks <- list()
for (i in starts) {
  label <- sub("^```\\{r ([^,}]+).*$", "\\1", lines[i])
  chunks[[label]] <- lines[(i + 1):(min(ends[ends > i]) - 1)]
}
run_chunk <- function(label) {
  cat("\n==== chunk:", label, "====", format(Sys.time()), "\n")
  eval(parse(text = chunks[[label]]), envir = globalenv())
}

run_chunk("libraries")
run_chunk("options")
RUN_REGRESSIONS    <- TRUE       # sample when the cache is missing
run_extra_analyses <- TRUE       # the fit chunk lives in Extra Analyses
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) bystudy_ranges <- args[1]   # one range only; the chunk defaults to both
run_chunk("model-helpers")       # load_or_fit()
run_chunk("read-data")
run_chunk("clean-trial")
run_chunk("boxlong")
run_chunk("flex-design")         # build_flex_dat(), util_cols(), dwell_cols()
run_chunk("extra-bystudy-fit")   # the long part: Narrow, then Wide
cat("\nDone", format(Sys.time()), "\n")
