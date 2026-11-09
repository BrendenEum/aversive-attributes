# Range normalization (Padoa-Schioppa 2009; Soltani, De Martino & Camerer 2012;
# Koszegi & Szeidl 2013). A box's subjective value is where its objective value
# sits in the range of its attribute type,
#   f(a) = (a - min) / (max - min),
# with min and max the smallest and largest value that type takes in the
# participant's range condition. The worst value of a type maps to 0 and the
# best to 1 in both conditions, so one objective point is worth 1/25 of the
# range in Narrow (P and N run from 5 to 30 in magnitude) but 1/45 in Wide
# (5 to 50). An option's subjective value is the sum of its two boxes.
#
# Shared pieces (trial sets, participants, logit choice rule) are in
# code/simulation/simulate_normalized_choices.R. Needs dplyr attached.
source(here::here("code/simulation/simulate_normalized_choices.R"))

# Subjective value of a vector of boxes, given each box's type and range
# condition. `design` comes from design_values().
range_normalize <- function(value, type, WN, design) {
  f <- numeric(length(value))
  for (i in seq_along(value)) {
    ctx  <- context_values(design, WN[i], type[i])
    f[i] <- (value[i] - min(ctx)) / (max(ctx) - min(ctx))
  }
  f
}

# Add each option's subjective value (the sum of its two boxes) to a trial set.
value_range_normalization <- function(trials, design) {
  trials$UL <- range_normalize(trials$VLT, trials$AttT, trials$WN, design) +
               range_normalize(trials$VLB, trials$AttB, trials$WN, design)
  trials$UR <- range_normalize(trials$VRT, trials$AttT, trials$WN, design) +
               range_normalize(trials$VRB, trials$AttB, trials$WN, design)
  trials
}

# Simulate every participant's choices on the trials of their condition.
simulate_range_normalization <- function(trials, participants, design, seed = 3) {
  simulate_logit_choices(value_range_normalization(trials, design), participants, seed)
}
