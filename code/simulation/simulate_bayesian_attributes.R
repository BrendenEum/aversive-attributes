# Bayesian attribute-level model of choice, after Fromer, Callaway, Griffiths &
# Shenhav (2025, Open Mind), for the 2-option x 2-attribute task.
#
# Every box (option x attribute) is an item whose value the participant
# estimates by sampling while the box is hovered. The estimate starts at a
# prior mean and moves toward what was seen; how far it moves depends on how
# long the box was looked at and on how noisy looking is.
#
# Two reference distributions, each either per attribute type or per range
# condition, say where the range enters:
#   units     the distribution whose spread sets the units of value and noise:
#             z = (value - mu_u) / sd_u, so a second of looking is equally
#             noisy in those units in both range conditions
#               "type":   the six values of the box's attribute type in that
#                         range condition (noise scales with the type's spread)
#               "cond":   all eighteen values of the range condition (same
#                         scaling for P, N and B)
#               "points": no rescaling (z = value; the control for noise scaling)
#   prior_ref the distribution the prior is calibrated to:
#               "type":   one prior per attribute type x range condition
#               "cond":   one prior per range condition, shared by P, N and B
#
# Prior. Normal(m0, prior_sd^2) in z units. Its mean sits a fraction q of the
# way from the prior reference's minimum to its maximum (q = 0.5 is the
# middle), and its SD is s0 times the prior reference's spread (s0 = 1 is
# calibrated). prior = "zero" puts the mean at an objective value of 0 instead
# (Fromer et al.'s ZERO PRIOR; with the symmetric spans of this design that is
# the middle, q = 0.5); prior = "flat" keeps the mean at the middle and makes
# s0 so large that the prior carries no weight in the estimate. The mean still
# matters under the flat prior: it is where the loss-aversion kink sits.
#
# Sampling. t_ij seconds of dwell on attribute j of option i give a sample mean
# xbar_ij ~ Normal(z_ij, sigma_x^2 / t_ij). The box's posterior mean is
# m_ij = m0 + pi(t_ij) * (xbar_ij - m0), where pi(t_ij) is the relative precision
# of the sampled information (the paper's notation; w in the code below),
#   pi(t_ij) = (t_ij / sigma_x^2) / (1 / sigma_p^2 + t_ij / sigma_x^2),
# with sigma_x = sigma and sigma_p = prior_sd in the code.
#
# Loss aversion. The posterior mean is valued relative to the prior mean: a box
# whose estimate falls below m0 counts lambda times its distance,
#   u(m) = m0 + (m - m0)            if m >= m0
#   u(m) = m0 + lambda * (m - m0)   if m <  m0,
# so with lambda > 1 attention repels more strongly from a below-reference box
# than it attracts toward an above-reference one. lambda = 1 is the plain model.
# The kink sits at the prior mean whatever the prior, so under the zero and
# flat priors it sits at 0.
# An option's estimate is the sum of its two boxes' u(m); the larger is chosen.
#
# Choice probability. Given the dwell times each m is normal, so each box's
# u(m) is a normal whose lower half is stretched by lambda. With lambda = 1
# nothing is stretched, V_left - V_right is normal, and P(choose left) =
# pnorm(Delta / sqrt(v)) with Delta and v its mean and variance (kinked_moments
# gives both in closed form for any lambda). With lambda != 1 the difference is
# not normal and that formula is only an approximation, so the probability is
# computed exactly on a grid of utility values instead (piece_choice_logprob
# below and code/simulation/exact_choice_prob.cpp). simulate_choices_direct()
# checks the probabilities against a box-by-box simulation.
#
# Parameters, collected in a list `par`: sigma (noise per second; one per
# subject, passed separately), s0, q and lambda (vectors named "Narrow" and
# "Wide", or one number for both). A "variant" is list(units, prior_ref,
# prior, lambda, ...); see bayesian_variants() in
# code/fitting/fit_bayesian_attributes.R.
#
# Needs dplyr attached (library(tidyverse)).

# ---------------------------------------------------------------------------
# Reference distributions and units
# ---------------------------------------------------------------------------

# One row per (range condition, attribute type) with the mean, SD, min and max
# of a reference distribution. Statistics are over the distinct design values
# (each value counted once), not over how often each value appeared.
reference_distributions <- function(boxlong, model = "type") {
  vals <- boxlong %>% distinct(WN, type, value)
  if (model == "type") {
    ref <- vals %>%
      group_by(WN, type) %>%
      summarise(mu = mean(value), sd = sd(value), min = min(value), max = max(value), .groups = "drop")
  } else if (model == "cond") {
    ref <- vals %>%
      group_by(WN) %>%
      summarise(mu = mean(value), sd = sd(value), min = min(value), max = max(value), .groups = "drop")
    # The same row for every attribute type.
    ref <- bind_rows(
      ref %>% mutate(type = "P"),
      ref %>% mutate(type = "N"),
      ref %>% mutate(type = "B")
    ) %>% mutate(type = factor(type, levels = levels(boxlong$type)))
  } else {
    stop("model must be 'type' or 'cond'")
  }
  ref %>% select(WN, type, mu, sd, min, max)
}

# Express every box in its units (z) and attach the prior reference's
# statistics (mu_p, sd_p, min_p, max_p) so the prior can be placed later.
box_units <- function(boxlong, variant) {
  ref_u <- reference_distributions(boxlong, if (variant$units == "points") "type" else variant$units) %>%
    transmute(WN, type, mu_u = mu, sd_u = sd)
  ref_p <- reference_distributions(boxlong, variant$prior_ref) %>%
    transmute(WN, type, mu_p = mu, sd_p = sd, min_p = min, max_p = max)
  boxes <- boxlong %>%
    left_join(ref_u, by = c("WN", "type")) %>%
    left_join(ref_p, by = c("WN", "type"))
  if (variant$units == "points") {        # no rescaling: z is the value in points
    boxes$mu_u <- 0
    boxes$sd_u <- 1
  }
  boxes$z        <- (boxes$value - boxes$mu_u) / boxes$sd_u
  boxes$cond_chr <- as.character(boxes$WN)   # for looking up per-condition parameters
  boxes
}

# Prior mean of every box, in z units. q is one number or one entry per box.
#   prior = "fitted": a fraction q of the way from the prior reference's
#                     minimum to its maximum
#   prior = "zero":   at an objective value of 0 (Fromer et al.'s ZERO PRIOR)
#   prior = "flat":   the middle of the prior reference (0 here); it only acts
#                     as the kink of the loss aversion and as the estimate of a
#                     box that was never seen
box_prior_mean <- function(boxes, q, prior = "fitted") {
  if (prior == "zero") {
    pos <- rep(0, nrow(boxes))
  } else {
    if (prior == "flat") q <- 0.5
    pos <- boxes$min_p + q * (boxes$max_p - boxes$min_p)   # in points
  }
  (pos - boxes$mu_u) / boxes$sd_u
}

# Prior SD of every box, in z units: s0 times the prior reference's spread.
box_prior_sd <- function(boxes, s0, prior = "fitted") {
  if (prior == "flat") s0 <- 1e3
  s0 * boxes$sd_p / boxes$sd_u
}

# Expand a per-condition parameter (one number, or a vector named "Narrow" /
# "Wide") to one entry per box. cond_chr is the box's condition as a string.
by_condition <- function(x, cond_chr) {
  if (length(x) == 1) rep(unname(x), length(cond_chr)) else unname(x[cond_chr])
}

# ---------------------------------------------------------------------------
# Bayesian updating of one box, in closed form
# ---------------------------------------------------------------------------

# Mean and variance of a box's posterior mean after `dwell` seconds of looking,
# given the prior Normal(m0, prior_sd^2) and noise sigma per second (all in z
# units). A box that was never seen keeps the prior mean, with no variance.
# sigma is one number, one value per box, or a matrix with one row per box and
# one column per candidate sigma (the fitting code evaluates every box at a
# whole grid of sigmas at once); the per-box inputs recycle down the columns.
bayes_box_moments <- function(z, dwell, m0, prior_sd, sigma) {
  seen  <- as.numeric(dwell > 0)
  tau_t <- dwell / sigma^2                             # precision contributed by the samples
  w     <- tau_t / (1 / prior_sd^2 + tau_t) * seen     # weight on the samples; 0 for a box never seen
  m_mean <- m0 + w * (z - m0)
  m_var  <- w^2 * sigma^2 / ifelse(seen == 1, dwell, 1)   # 0 for a box never seen, since w = 0 there
  list(mean = m_mean, var = m_var)
}

# ---------------------------------------------------------------------------
# Pieces: the data structure the fitting and simulation code works on
# ---------------------------------------------------------------------------

# A "piece" is a set of boxes with a local integer trial_id (1..n) plus the
# matching trials table (one row per trial_id, in that order). prepare_boxes()
# builds, for one variant:
#   $all        one piece with every box; its trial ids are the global ids
#   $cond[[c]]  one piece per range condition ("Narrow", "Wide")
#   $subj[[s]]  one piece per subject
# Every trials table carries global_id (the trial's id in $all) so choices and
# train/test masks defined on $all$trials can be mapped to any piece.
# Trials with a box whose dwell time is NA (never hovered / not recorded) are
# dropped, as the paper's regressions drop them; $n_dropped says how many.
prepare_boxes <- function(boxlong, variant) {
  boxes <- box_units(boxlong, variant)
  key   <- paste(boxes$Subject, boxes$Trial)
  drop  <- key %in% unique(key[is.na(boxes$dwell)])
  n_dropped <- length(unique(key[drop]))
  boxes <- boxes[!drop, ]
  key   <- key[!drop]
  boxes$global_id <- as.integer(factor(key, levels = unique(key)))

  make_piece <- function(bx) {
    bx$trial_id <- as.integer(factor(bx$global_id, levels = unique(bx$global_id)))
    trials <- bx[!duplicated(bx$global_id), c("Study", "Subject", "Trial", "WN", "ChooseL", "global_id")]
    trials$trial_id <- seq_len(nrow(trials))
    rownames(trials) <- NULL
    list(boxes = bx, trials = trials)
  }

  all  <- make_piece(boxes)
  cond <- list()
  for (cc in levels(boxes$WN)) cond[[cc]] <- make_piece(boxes[boxes$WN == cc, ])
  subj <- list()
  by_subject <- split(boxes, boxes$Subject)
  for (s in names(by_subject)) subj[[s]] <- make_piece(by_subject[[s]])

  list(all = all, cond = cond, subj = subj, variant = variant,
       subjects = all$trials %>% distinct(Study, Subject, WN), n_dropped = n_dropped)
}

# Overwrite the observed choices in every piece of a prep with new ones (e.g.
# simulated), given as a 0/1 vector aligned with prep$all$trials (global_id order).
set_choices <- function(prep, choices_by_global_id) {
  prep$all$trials$ChooseL <- choices_by_global_id[prep$all$trials$global_id]
  for (cc in names(prep$cond)) prep$cond[[cc]]$trials$ChooseL <- choices_by_global_id[prep$cond[[cc]]$trials$global_id]
  for (s  in names(prep$subj)) prep$subj[[s]]$trials$ChooseL  <- choices_by_global_id[prep$subj[[s]]$trials$global_id]
  prep
}

# Mean and variance of u(m) for every box, where m ~ Normal(m_mean, m_var) is
# the posterior mean and u applies the loss-aversion kink at the prior mean m0.
# Partial moments of X = m - m0: with a = mean / sd,
#   E[X+]   = mean * Phi(a) + sd * phi(a)          E[X-]   = mean - E[X+]
#   E[X+^2] = (mean^2 + sd^2) Phi(a) + mean sd phi(a)
#   E[X-^2] = (mean^2 + sd^2) (1 - Phi(a)) - mean sd phi(a)
# and u(m) - m0 = X+ + lambda * X-, with X+ X- = 0.
kinked_moments <- function(m_mean, m_var, m0, lambda) {
  mu <- m_mean - m0
  sd <- sqrt(pmax(m_var, 0))
  a  <- ifelse(sd > 0, mu / sd, ifelse(mu >= 0, Inf, -Inf))
  Phi <- pnorm(a)
  phi <- dnorm(a)
  EXp  <- mu * Phi + sd * phi
  EXm  <- mu - EXp
  EX2p <- (mu^2 + sd^2) * Phi + mu * sd * phi
  EX2m <- (mu^2 + sd^2) * (1 - Phi) - mu * sd * phi
  u_mean <- m0 + EXp + lambda * EXm
  u_var  <- pmax(EX2p + lambda^2 * EX2m - (EXp + lambda * EXm)^2, 0)
  list(mean = u_mean, var = u_var)
}

# Mean and variance of V_left - V_right for every trial of a piece. sigma is a
# per-box vector or one number (vectors out), or a boxes x candidates matrix
# (trials x candidates matrices out); par has s0, q and lambda (per condition).
piece_moments <- function(piece, variant, sigma, par) {
  b   <- piece$boxes
  m0  <- box_prior_mean(b, by_condition(par$q, b$cond_chr), variant$prior)
  psd <- box_prior_sd(b, by_condition(par$s0, b$cond_chr), variant$prior)
  mom <- bayes_box_moments(b$z, b$dwell, m0, psd, sigma)
  lam <- by_condition(par$lambda, b$cond_chr)
  if (any(lam != 1)) mom <- kinked_moments(mom$mean, mom$var, m0, lam)
  delta <- rowsum(b$side_sign * mom$mean, b$trial_id)
  var   <- rowsum(mom$var, b$trial_id)
  if (is.matrix(sigma)) list(delta = delta, var = var) else list(delta = as.vector(delta), var = as.vector(var))
}

# ---------------------------------------------------------------------------
# Choice probabilities
# ---------------------------------------------------------------------------

# With lambda = 1, V_left - V_right is normal and P(choose left) is the probit
# pnorm(Delta / sqrt(v)). With lambda != 1 each box's u(m) is a normal whose
# lower half is stretched, the difference is not normal, and the probit is only
# an approximation (off by up to 0.006 at the fitted lambda of 1.065, and by up
# to 0.08 at lambda = 2.5). The probability is then computed exactly, in C++
# (exact_choice_prob.cpp): every box's utility is laid on a grid of values,
# each option's utility is the sum of its boxes, and the left option is
# compared with the right cell by cell. The grid's error shrinks with the
# square of its spacing, so the calculation is run at spacings h, h/2 and h/4
# and extrapolated to zero spacing (Richardson extrapolation), which removes
# the h^2 and h^4 error terms. Checked against the probit at lambda = 1
# (agreement to 2e-9), against 1e7-draw simulations at lambda = 1.065 and 2.5
# (agreement within the simulation error, 5e-5), and against a grid eight
# times finer (relative differences 2e-7 for ordinary probabilities, 1e-5
# around 1e-8; the full data's log-likelihood moves by ~3e-5 at the fitted
# parameters, the number the supplement reports).
#
# The probit stays exact when no box of the trial can land on its stretched
# side, i.e. every box is more than KINK_SD SDs from the reference point (a
# probability below 1e-16 of landing there); those trials keep it, which is
# cheaper. Boxes that were never seen (kappa = 0) sit exactly at the reference
# point and add nothing random; the C++ leaves them out of the sum.
GRID_FACTOR <- 0.5           # grid spacing h = GRID_FACTOR * the trial's smallest box SD (kappa)
GRID_NSD    <- 7.5           # each box's grid covers this many SDs on either side
KINK_SD     <- 8.3           # pnorm(-8.3) < 1e-16: beyond this the kink carries no probability
Rcpp::sourceCpp(here::here("code/simulation/exact_choice_prob.cpp"))

# Extrapolation of the three grid levels (spacings h, h/2, h/4) to zero
# spacing: the levels' errors are a h^2 + b h^4 + ..., and these weights cancel
# a and b.
extrapolate_levels <- function(P) (P[, 1] - 20 * P[, 2] + 64 * P[, 3]) / 45

# Log probabilities of choosing left and of choosing right for every trial of
# a piece, as list(left, right): trials x candidates matrices when sigma is a
# boxes x candidates matrix, vectors when sigma is a per-box vector or one
# number. par has s0, q and lambda (per condition). Both sides are computed
# directly (the right one is not 1 minus the left), so a small probability
# keeps its relative accuracy and its log is sound. exact = FALSE returns the
# probit for any lambda, only for showing the size of the approximation.
piece_choice_logprob <- function(piece, variant, sigma, par, exact = TRUE) {
  b   <- piece$boxes
  vec <- !is.matrix(sigma)
  if (vec) sigma <- matrix(sigma, nrow(b), 1)
  m0  <- box_prior_mean(b, by_condition(par$q, b$cond_chr), variant$prior)
  psd <- box_prior_sd(b, by_condition(par$s0, b$cond_chr), variant$prior)
  lam <- by_condition(par$lambda, b$cond_chr)
  mom <- bayes_box_moments(b$z, b$dwell, m0, psd, sigma)       # boxes x candidates
  kinked <- any(lam != 1)
  km    <- if (kinked) kinked_moments(mom$mean, mom$var, m0, lam) else mom
  delta <- rowsum(b$side_sign * km$mean, b$trial_id)            # trials x candidates
  v     <- rowsum(km$var, b$trial_id)
  z     <- ifelse(v > 0, delta / sqrt(v), 0)
  left  <- pnorm(z, log.p = TRUE)       # the log of the normal CDF directly, so confident trials don't underflow
  right <- pnorm(-z, log.p = TRUE)
  if (kinked && exact) {
    nt <- nrow(delta)
    stopifnot(nrow(b) == 4 * nt, all(tabulate(b$trial_id, nt) == 4))
    eta_b <- mom$mean - m0                                       # each box's estimate relative to the reference point: mean
    kap_b <- sqrt(mom$var)                                       # and SD
    # Trials x candidates where some box can land on its stretched side (and
    # lambda != 1): the probit is an approximation there and the grid takes over.
    near  <- kap_b > 0 & abs(eta_b) <= KINK_SD * kap_b
    lam_t <- lam[match(seq_len(nt), b$trial_id)]                 # each trial's lambda (its boxes share it)
    need  <- which(rowsum(near + 0, b$trial_id) > 0 & lam_t != 1)   # trials x candidates (near + 0 keeps the matrix)
    if (length(need) > 0) {
      rows <- matrix(order(b$trial_id), nrow = 4)                # column t: the four boxes of trial t
      col4 <- function(x) matrix(sapply(1:4, function(p) as.vector(x[rows[p, ], , drop = FALSE])), ncol = 4)
      ng   <- ncol(sigma)
      P <- exact_choice_prob(col4(eta_b)[need, , drop = FALSE], col4(kap_b)[need, , drop = FALSE],
                             col4(matrix(b$side_sign, nrow(b), ng))[need, , drop = FALSE],
                             col4(matrix(lam, nrow(b), ng))[need, , drop = FALSE],
                             GRID_FACTOR, 3L, GRID_NSD)
      left[need]  <- log(pmax(extrapolate_levels(P[, 1:3, drop = FALSE]), 0))
      right[need] <- log(pmax(extrapolate_levels(P[, 4:6, drop = FALSE]), 0))
    }
  }
  if (vec) list(left = as.vector(left), right = as.vector(right)) else list(left = left, right = right)
}

# P(choose left) for every trial of a piece (a trials x candidates matrix when
# sigma is a matrix, a vector otherwise).
piece_choice_prob <- function(piece, variant, sigma, par, exact = TRUE) {
  exp(piece_choice_logprob(piece, variant, sigma, par, exact)$left)
}

# ---------------------------------------------------------------------------
# Simulating choices
# ---------------------------------------------------------------------------

# P(choose left) for every trial of the full data, as a data frame aligned with
# prep$all$trials. sigma_by_subject is a vector named by Subject.
choice_probabilities <- function(prep, sigma_by_subject, par) {
  piece <- prep$all
  sigma <- unname(sigma_by_subject[as.character(piece$boxes$Subject)])
  out <- piece$trials
  out$p_left <- piece_choice_prob(piece, prep$variant, sigma, par)
  out
}

# Draw `reps` simulated choices per trial from the choice probabilities. Each
# rep gets its own trial id (Trial_rep) so reps can be stacked and treated as
# separate trials by the regression code.
simulate_choices_from_prob <- function(trials, reps = 10, seed = 1) {
  set.seed(seed)
  out <- list()
  for (r in seq_len(reps)) {
    d <- trials
    d$rep         <- r
    d$Trial_rep   <- d$Trial + 1e6 * r
    d$ChooseL_sim <- as.integer(runif(nrow(d)) < d$p_left)
    out[[r]] <- d
  }
  bind_rows(out)
}

# The same model simulated box by box (draw every box's sample mean, update,
# apply the kink, sum, compare). Only used to check the closed form.
simulate_choices_direct <- function(prep, sigma_by_subject, par, seed = 1) {
  set.seed(seed)
  piece <- prep$all
  b     <- piece$boxes
  v     <- prep$variant
  sigma <- unname(sigma_by_subject[as.character(b$Subject)])
  m0    <- box_prior_mean(b, by_condition(par$q, b$cond_chr), v$prior)
  psd   <- box_prior_sd(b, by_condition(par$s0, b$cond_chr), v$prior)
  lam   <- by_condition(par$lambda, b$cond_chr)
  tau_t <- b$dwell / sigma^2
  w     <- tau_t / (1 / psd^2 + tau_t)
  w[b$dwell <= 0] <- 0
  sample_sd <- ifelse(b$dwell > 0, sigma / sqrt(pmax(b$dwell, 1e-9)), 0)
  xbar  <- b$z + rnorm(nrow(b), 0, sample_sd)
  m     <- m0 + w * (xbar - m0)
  u     <- m0 + (m - m0) * ifelse(m >= m0, 1, lam)
  V     <- as.vector(rowsum(b$side_sign * u, b$trial_id))
  out <- piece$trials
  out$ChooseL_sim <- as.integer(V > 0)
  out
}
