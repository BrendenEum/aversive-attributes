# Sample the two attention-weight models (Narrow, Wide; each pooled over the
# studies) outside of RStudio, so that a closed session can't kill a run that
# takes many hours. From the project root:
#
#   caffeinate -i nohup Rscript code/fitting/sample_pooled_attention_weights.R > results/logs/sampling.log 2>&1 &
#
# With a range as the argument it samples only that fit, so the two can run at
# the same time (four chains each, eight cores):
#
#   caffeinate -i nohup Rscript code/fitting/sample_pooled_attention_weights.R Narrow > results/logs/sampling_narrow.log 2>&1 &
#   caffeinate -i nohup Rscript code/fitting/sample_pooled_attention_weights.R Wide   > results/logs/sampling_wide.log   2>&1 &
#
# It evaluates the notebook's own chunks (libraries, data, design matrices,
# fits), so the caches it leaves in notebooks/brm_fits/ are exactly what
# Main_Analysis.Rmd loads afterwards. brm() skips any fit whose cache already
# exists, so the script can be re-run after an interruption.

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
RUN_REGRESSIONS <- TRUE          # sample when the cache is missing
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) fit_ranges <- args[1]   # one range only; the chunk defaults to both
run_chunk("model-helpers")       # load_or_fit()
run_chunk("read-data")
run_chunk("clean-trial")
run_chunk("boxlong")
run_chunk("flex-design")
run_chunk("flex-fit-by-range")   # the long part: Narrow, then Wide
cat("\nDone", format(Sys.time()), "\n")
