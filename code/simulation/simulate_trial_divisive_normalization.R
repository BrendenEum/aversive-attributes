# Divisive normalization within the trial (Louie, Khaw & Glimcher 2013; Landry
# & Webb 2021, "pairwise normalization"). Each attribute is normalized by the
# two values it takes on the current trial, one per option:
#   f(a_L) = a_L / (sigma_h + omega * (|a_L| + |a_R|)),   and likewise for a_R,
# so a box's subjective value is its share of the attribute's total magnitude
# on that trial. The range condition never enters directly: the Wide condition
# is less sensitive to objective differences only because its trials are built
# from larger values, which make for larger denominators. An option's
# subjective value is the sum of its two boxes.
#
# Parameters as in the context version: sigma_h = 5 (the smallest magnitude in
# the design) and omega = 1. With sigma_h = 0 the smallest mixed-sign pair
# (5 vs -5) would be exactly as discriminable as the largest (15 vs -15);
# sigma_h > 0 softens that.
#
# Shared pieces (trial sets, participants, logit choice rule) are in
# code/simulation/simulate_normalized_choices.R. Needs dplyr attached.
source(here::here("code/simulation/simulate_normalized_choices.R"))

# Subjective values of the left and right box of one attribute, given the two
# objective values. Vectorized over trials.
trial_divisive_normalize <- function(vL, vR, sigma_h = 5, omega = 1) {
  denom <- sigma_h + omega * (abs(vL) + abs(vR))
  list(L = vL / denom, R = vR / denom)
}

# Add each option's subjective value (the sum of its two boxes) to a trial set.
value_trial_divisive_normalization <- function(trials, sigma_h = 5, omega = 1) {
  top    <- trial_divisive_normalize(trials$VLT, trials$VRT, sigma_h, omega)
  bottom <- trial_divisive_normalize(trials$VLB, trials$VRB, sigma_h, omega)
  trials$UL <- top$L + bottom$L
  trials$UR <- top$R + bottom$R
  trials
}

# Simulate every participant's choices on the trials of their condition.
# `design` is not used by this model; it is accepted so that every model can be
# called the same way.
simulate_trial_divisive_normalization <- function(trials, participants, design = NULL, seed = 3, sigma_h = 5, omega = 1) {
  simulate_logit_choices(value_trial_divisive_normalization(trials, sigma_h, omega), participants, seed)
}
