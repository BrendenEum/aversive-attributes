# Divisive normalization by the context (Louie, Khaw & Glimcher 2013; Khaw,
# Glimcher & Louie 2017; Zimmermann, Glimcher & Louie 2018). A box's objective
# value is divided by the typical magnitude of the values its attribute type
# takes in the participant's range condition,
#   f(a) = a / (sigma_h + omega * mean(|context|)),
# where sigma_h is a semi-saturation constant and omega weights the context.
# The denominator is a property of the condition, so within a condition the
# transform is a straight line, but its slope is shallower in Wide, where the
# values are larger on average (mean |a| is 17.5 in Narrow and 25.8 in Wide
# for P and N). An option's subjective value is the sum of its two boxes.
#
# Parameters: sigma_h = 5 (the smallest magnitude in the design) and
# omega = 1. The Narrow > Wide prediction holds for any sigma_h >= 0; it only
# fades toward equality as sigma_h grows and swamps the context.
#
# Shared pieces (trial sets, participants, logit choice rule) are in
# code/simulation/simulate_normalized_choices.R. Needs dplyr attached.
source(here::here("code/simulation/simulate_normalized_choices.R"))

# Subjective value of a vector of boxes, given each box's type and range
# condition. `design` comes from design_values().
context_divisive_normalize <- function(value, type, WN, design, sigma_h = 5, omega = 1) {
  f <- numeric(length(value))
  for (i in seq_along(value)) {
    ctx  <- context_values(design, WN[i], type[i])
    f[i] <- value[i] / (sigma_h + omega * mean(abs(ctx)))
  }
  f
}

# Add each option's subjective value (the sum of its two boxes) to a trial set.
value_context_divisive_normalization <- function(trials, design, sigma_h = 5, omega = 1) {
  trials$UL <- context_divisive_normalize(trials$VLT, trials$AttT, trials$WN, design, sigma_h, omega) +
               context_divisive_normalize(trials$VLB, trials$AttB, trials$WN, design, sigma_h, omega)
  trials$UR <- context_divisive_normalize(trials$VRT, trials$AttT, trials$WN, design, sigma_h, omega) +
               context_divisive_normalize(trials$VRB, trials$AttB, trials$WN, design, sigma_h, omega)
  trials
}

# Simulate every participant's choices on the trials of their condition.
simulate_context_divisive_normalization <- function(trials, participants, design, seed = 3, sigma_h = 5, omega = 1) {
  simulate_logit_choices(value_context_divisive_normalization(trials, design, sigma_h, omega), participants, seed)
}
