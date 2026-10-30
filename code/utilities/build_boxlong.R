# Build the trial-level (`pas`) and box-level (`boxlong`) data frames from the
# merged data in data/processed/main_data.csv. This is the same cleaning that
# Main_Analysis.Rmd does in its `clean-trial` and `boxlong` chunks, moved into
# two functions so that scripts under code/ (fitting, parameter recovery) and
# the exploratory notebooks build exactly the same data. getAllUtilities.R
# sources this file automatically. Needs dplyr attached (library(tidyverse)).

# One row per trial. Narrow is the reference range condition; ChooseL and
# Correct become 0/1; values are also kept in units of 10 points because the
# regressions use that scale (the paper's SV_* convention).
build_pas <- function(pas_raw) {
  pas_raw %>%
    mutate(
      WN      = factor(WN,  levels = c("Narrow", "Wide")),
      BPN     = factor(BPN, levels = c("PN", "PB", "NB")),
      ChooseL = as.integer(ChooseL),
      Correct = as.integer(Correct),
      Difficulty   = abs(Diff),          # best minus worst option value
      Diff10       = Diff / 10,
      Difficulty10 = Difficulty / 10
    )
}

# One row per box (left/right option x top/bottom attribute), with the box's
# attribute type, value, and dwell time in seconds. Left boxes get a side sign
# of +1 (looking there should push choice toward "left"), right boxes -1.
build_boxlong <- function(pas) {
  mk_box <- function(df, box, side_sign, type_col, val_col, dwell_col) {
    df %>%
      transmute(
        Study, Subject, Trial, WN, ChooseL,
        box       = box,
        side_sign = side_sign,
        type      = .data[[type_col]],
        value     = .data[[val_col]],
        dwell     = .data[[dwell_col]]
      )
  }
  bind_rows(
    mk_box(pas, "LT",  1, "AttT", "VLT", "DwellLT"),
    mk_box(pas, "LB",  1, "AttB", "VLB", "DwellLB"),
    mk_box(pas, "RT", -1, "AttT", "VRT", "DwellRT"),
    mk_box(pas, "RB", -1, "AttB", "VRB", "DwellRB")
  ) %>%
    mutate(
      type = factor(type, levels = c("P", "N", "B")),
      # Signed dwell: seconds, positive if it should push choice toward this box's side.
      signed_dwell = side_sign * dwell,
      # Signed value, in units of 10 points.
      signed_value = side_sign * value / 10,
      # Cell label for (type x value), e.g. "N_neg10", "P_pos5". Used to name
      # the attention-weight coefficients.
      cell = paste0(type, "_", if_else(value < 0, "neg", "pos"), abs(value))
    )
}
