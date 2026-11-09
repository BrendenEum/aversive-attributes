# Rank normalization, after decision by sampling (Stewart, Chater & Brown 2006;
# Walasek & Stewart 2015). A box's subjective value is the rank of its
# objective value among the values its attribute type takes in the
# participant's range condition, rescaled to run from 0 (the worst value) to 1
# (the best):
#   f(a) = (rank(a) - 1) / (n - 1),
# with n = 6 values per type. This is the empirical distribution function of
# the context (Bhui & Gershman 2018). Both conditions have six values per
# type, so one rank step is worth 5 objective points in Narrow but 9 on
# average in Wide, and subjective value grows more slowly per objective point
# in Wide. Rank and range normalization coincide when the values are equally
# spaced (Narrow) and differ when they are not (Wide). An option's subjective
# value is the sum of its two boxes.
#
# Shared pieces (trial sets, participants, logit choice rule) are in
# code/simulation/simulate_normalized_choices.R. Needs dplyr attached.
source(here::here("code/simulation/simulate_normalized_choices.R"))

# Subjective value of a vector of boxes, given each box's type and range
# condition. `design` comes from design_values().
rank_normalize <- function(value, type, WN, design) {
  f <- numeric(length(value))
  for (i in seq_along(value)) {
    ctx  <- context_values(design, WN[i], type[i])
    f[i] <- sum(ctx < value[i]) / (length(ctx) - 1)      # (rank - 1) / (n - 1)
  }
  f
}

# Add each option's subjective value (the sum of its two boxes) to a trial set.
value_rank_normalization <- function(trials, design) {
  trials$UL <- rank_normalize(trials$VLT, trials$AttT, trials$WN, design) +
               rank_normalize(trials$VLB, trials$AttB, trials$WN, design)
  trials$UR <- rank_normalize(trials$VRT, trials$AttT, trials$WN, design) +
               rank_normalize(trials$VRB, trials$AttB, trials$WN, design)
  trials
}

# Simulate every participant's choices on the trials of their condition.
simulate_rank_normalization <- function(trials, participants, design, seed = 3) {
  simulate_logit_choices(value_rank_normalization(trials, design), participants, seed)
}
