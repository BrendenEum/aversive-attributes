# Shared pieces of the range-sensitivity simulations. Each of the four scripts
#   code/simulation/simulate_range_normalization.R
#   code/simulation/simulate_context_divisive_normalization.R
#   code/simulation/simulate_trial_divisive_normalization.R
#   code/simulation/simulate_rank_normalization.R
# defines one model of value normalization. The pieces here are the same for
# every model: the two hypothetical studies, the simulated participants, and the
# logit choice rule.
#
# The question these simulations answer: if subjective value is normalized to
# the values a participant has seen, is choice less sensitive to a one-point
# difference in objective value when the range is Wide than when it is Narrow?
# Every model maps a box's objective value to a subjective value through the
# context of its attribute type (the six values that type can take in the
# participant's range condition, or, for the trial version of divisive
# normalization, the two values on the current trial). An option's subjective
# value is the sum of its two boxes, and the probability of choosing left is a
# logistic function of the left minus right subjective value.
#
# Needs dplyr attached (library(tidyverse)) and `pas`, the trial-level data
# frame of notebooks/Main_Analysis.Rmd (one row per trial, with the type and
# value of each box: AttT, AttB, VLT, VRT, VLB, VRB).

# ---------------------------------------------------------------------------
# Context
# ---------------------------------------------------------------------------

# The values each attribute type (P, N, B) can take in each range condition,
# read off the real trials. Returns one row per (WN, type, value); six values
# per type and condition in this design.
design_values <- function(pas) {
  boxes <- rbind(
    data.frame(WN = pas$WN, type = pas$AttT, value = pas$VLT),
    data.frame(WN = pas$WN, type = pas$AttT, value = pas$VRT),
    data.frame(WN = pas$WN, type = pas$AttB, value = pas$VLB),
    data.frame(WN = pas$WN, type = pas$AttB, value = pas$VRB)
  )
  boxes <- unique(boxes)
  boxes <- boxes[order(boxes$WN, boxes$type, boxes$value), ]
  rownames(boxes) <- NULL
  boxes
}

# The values in one box's context: the design values of its type in its range
# condition.
context_values <- function(design, WN, type) {
  design$value[as.character(design$WN) == as.character(WN) & design$type == type]
}

# ---------------------------------------------------------------------------
# Hypothetical studies
# ---------------------------------------------------------------------------

# Two hypothetical studies, one per range condition, each made of 60 real
# trials: 20 from each block type (PN, PB, NB), drawn at random from all the
# trials of that condition in both real studies. Every simulated participant in
# a condition sees the same 60 trials. Each trial keeps the type and value of
# its four boxes (left/right x top/bottom) and the objective value difference
# between the options (Diff = left minus right, in points).
make_trial_sets <- function(pas, n_per_block = 20, seed = 1) {
  set.seed(seed)
  keep   <- c("WN", "BPN", "AttT", "AttB", "VLT", "VRT", "VLB", "VRB")
  trials <- data.frame()
  for (wn in c("Narrow", "Wide")) {
    for (bpn in c("PN", "PB", "NB")) {
      pool   <- pas[pas$WN == wn & pas$BPN == bpn, keep]
      picked <- pool[sample(nrow(pool), n_per_block), ]
      trials <- rbind(trials, picked)
    }
  }
  trials$Trial <- rep(1:(3 * n_per_block), times = 2)       # 1..60 within each condition
  trials$Diff  <- (trials$VLT + trials$VLB) - (trials$VRT + trials$VRB)
  rownames(trials) <- NULL
  trials
}

# ---------------------------------------------------------------------------
# Participants
# ---------------------------------------------------------------------------

# Simulated participants, n_per_condition in each range condition. Each has an
# inverse temperature beta (how sharply choice follows subjective value) and
# nothing else: no side bias, so a participant is indifferent between two
# options of equal subjective value. beta is drawn uniformly between beta_min
# and beta_max (1 to 5, mean 3). The range keeps every participant away from
# random responding (beta near 0) and from perfect choice on the normalized
# scales below, where one option's subjective value spans about 0 to 2.
# Nothing about a participant depends on their range condition, so any
# Narrow / Wide difference in the simulated choices comes from the
# normalization alone.
draw_participants <- function(n_per_condition = 260, beta_min = 1, beta_max = 5, seed = 2) {
  set.seed(seed)
  n <- 2 * n_per_condition
  data.frame(
    Subject = sprintf("Sim_%03d", 1:n),
    WN      = factor(rep(c("Narrow", "Wide"), each = n_per_condition), levels = c("Narrow", "Wide")),
    beta    = runif(n, min = beta_min, max = beta_max)
  )
}

# ---------------------------------------------------------------------------
# Choice rule
# ---------------------------------------------------------------------------

# Logit choice rule. `valued` is a trial set (from make_trial_sets) with the
# subjective value of each option, UL and UR, already filled in by a model.
# Every participant goes through the 60 trials of their condition and chooses
# left with probability
#   Pr(Left) = logistic(beta * (UL - UR)).
# The same seed for every model means the models are compared on the same
# random draws. Returns one row per participant x trial, with the regression
# variables of the notebook (ChooseL, Diff10 = Diff / 10).
simulate_logit_choices <- function(valued, participants, seed = 3) {
  set.seed(seed)
  sims <- list()
  for (i in 1:nrow(participants)) {
    p <- participants[i, ]
    t <- valued[valued$WN == p$WN, ]
    t$Subject <- p$Subject
    t$pLeft   <- plogis(p$beta * (t$UL - t$UR))
    t$ChooseL <- rbinom(nrow(t), size = 1, prob = t$pLeft)
    sims[[i]] <- t
  }
  out <- bind_rows(sims)
  out$Diff10 <- out$Diff / 10
  out
}
