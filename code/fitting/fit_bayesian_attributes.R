# Maximum-likelihood fitting of the Bayesian attribute-level model (the model
# itself is in code/simulation/simulate_bayesian_attributes.R; source that
# first). Choice probabilities come from piece_choice_logprob() there: the
# probit where it is exact (lambda = 1, so stages 1 and 2 below) and the exact
# grid calculation otherwise (stage 3 and the held-out predictions).
#
# Parameters
#   sigma   noise per second, one per subject. Not a fitted parameter but a
#           random effect: each subject's sigma is a draw from the population
#           distribution of their range condition (see Error model below)
#   mu, tau mean and SD of log sigma in the population, one pair per range
#           condition, fitted
#   q       prior position (fraction of the prior reference's range), one per range
#           condition; fixed at 0.5 for the zero and flat priors (see below)
#   lambda  loss aversion on the posterior relative to the prior mean: one value
#           shared by both range conditions (variant$lambda = "shared"), one per
#           condition ("cond"), or fixed at 1 ("fixed1", the no-loss-aversion checks)
# The prior SD is fixed at the spread of the prior reference (s0 = 1 in the
# model's parameterisation, 1e3 for the flat prior; an earlier round fitted s0
# freely and found 1.07 / 1.17). q and lambda are pooled across subjects: 30-60
# trials per person cannot pin down attention effects of ~0.1-0.3 log-odds per
# second of dwell.
#
# Error model. log sigma ~ Normal(mu, tau^2) in the subject's range condition,
# and a subject's contribution to the likelihood is the probability of their
# choices averaged over that distribution (the marginal likelihood, as in a
# mixed-effects regression). The average is a sum over a fine grid of log sigma
# (SIGMA_GRID) with weights proportional to the normal density at the grid
# points, so the population is a normal on log sigma truncated to the grid's
# range. This replaces an earlier penalised fit of one sigma per subject: plain
# maximum likelihood sends the sigma of a subject with 30 correct trials to 0,
# and a held-out error then costs 25+ log-likelihood units; averaging over the
# population puts little weight on those sigmas, and how much a subject's own
# trials pull away from the population is decided by the data (tau), not by a
# hand-set penalty. For reporting, a subject's sigma is the posterior mean of
# sigma given their trials and the fitted population; predictions average the
# choice probability over that posterior (posterior_choice_probabilities).
#
# One algorithm fits every variant; the variants differ only in which group
# parameters are free. Given q and lambda, the subjects are independent and
# each condition's (mu, tau) enters only that condition's subjects, so the fit
# works on the profile likelihood (fit_profile): at every candidate (q,
# lambda), each subject's log-likelihood is evaluated on the sigma grid once
# (grid_pass; the group parameters then only enter through the weights), and
# each condition's (mu, tau) is found by its own Nelder-Mead on those weighted
# sums (fit_population). Around that, three stages: (1) for a grid of q with
# lambda = 1 the NLL is recorded; (2) q is refined per condition from the best
# grid point by a bounded 1-D search (together with a per-condition lambda, by
# Nelder-Mead, when lambda = "cond"); (3) for a shared lambda, q in both
# conditions and lambda are refined jointly by Nelder-Mead. The population
# step runs from scratch at every candidate. When q is fixed (zero and flat
# priors) stage 1 is a single pass, stage 2 has nothing to refine and stage 3
# is a bounded 1-D search over lambda. The cross-validation (loo_by_trial_index:
# every subject's trials cut into blocks of TRIALS_PER_FOLD, fold k holding out
# block k) fits every fold from scratch by the same three stages on the
# training trials, population distribution included; nothing is carried over
# from the full fit.
#
# Variants: prior = "zero" fixes the prior position at an objective 0 (q = 0.5,
# the middle of the symmetric span); prior = "flat" keeps that position and
# sets s0 = 1e3, so the prior carries no weight in the estimate but the kink of
# the loss aversion still sits at 0. Both fit the population and lambda.
#
# Group parameters travel in a list `par` with s0, q, lambda, each a vector
# named by condition. Train/test masks are logical vectors over
# prep$all$trials (global_id order).

# The grid of log sigma the average over sigma runs on: 60 points from 0.005
# to 500 points per second (step 0.195). With the earlier probit likelihood,
# 60, 100, 200 and 400 points gave marginal NLLs within 1e-4 of each other and
# (mu, tau) equal to the fourth decimal; the supplement reports the same check
# (60 vs 120 points) at the fitted parameters with the exact likelihood.
SIGMA_GRID_LOG <- seq(log(0.005), log(500), length.out = 60)
SIGMA_GRID     <- exp(SIGMA_GRID_LOG)
MU_BOUNDS      <- range(SIGMA_GRID_LOG)
TAU_BOUNDS     <- c(0.2, 5)       # SD of log sigma across subjects; the floor is about one grid step (fitted values ~0.5, two and a half steps)
Q_BOUNDS       <- c(-1, 2)        # q may sit outside the reference range
LAMBDA_BOUNDS  <- c(0.25, 8)      # loss aversion (log scale); below 1 = gain seeking
INNER_CONTROL  <- list(maxit = 500, reltol = 1e-10)   # the population step: far tighter than the outer Nelder-Mead (1e-7)

# Nelder-Mead has no bounds, so the group parameters are optimised on an
# unconstrained scale and mapped into their bounds with a logistic: without
# this, a flat likelihood (e.g. a weak prior, where q no longer matters) lets
# the simplex drift to absurd values. mu and tau go through the same map.
to_theta   <- function(x, lo, hi) { p <- (x - lo) / (hi - lo); log(p / (1 - p)) }
from_theta <- function(th, lo, hi) lo + (hi - lo) / (1 + exp(-th))
th_q   <- function(q)   to_theta(q, Q_BOUNDS[1], Q_BOUNDS[2])
q_th   <- function(th)  from_theta(th, Q_BOUNDS[1], Q_BOUNDS[2])
th_lam <- function(lam) to_theta(log(lam), log(LAMBDA_BOUNDS[1]), log(LAMBDA_BOUNDS[2]))
lam_th <- function(th)  exp(from_theta(th, log(LAMBDA_BOUNDS[1]), log(LAMBDA_BOUNDS[2])))
th_mu  <- function(mu)  to_theta(mu, MU_BOUNDS[1], MU_BOUNDS[2])
mu_th  <- function(th)  from_theta(th, MU_BOUNDS[1], MU_BOUNDS[2])
th_tau <- function(tau) to_theta(log(tau), log(TAU_BOUNDS[1]), log(TAU_BOUNDS[2]))
tau_th <- function(th)  exp(from_theta(th, log(TAU_BOUNDS[1]), log(TAU_BOUNDS[2])))
clamp  <- function(x, lo, hi) pmin(pmax(x, lo + 1e-6 * (hi - lo)), hi - 1e-6 * (hi - lo))

# Negative log-likelihood of a piece's observed choices at given sigmas (a
# per-box vector or one number; only the trials where `use` is TRUE, NULL =
# all). Kept for checks; the fit works on the grid version below.
piece_nll <- function(piece, variant, sigma, par, use = NULL) {
  lp <- piece_choice_logprob(piece, variant, sigma, par)
  ll <- ifelse(piece$trials$ChooseL == 1, lp$left, lp$right)   # log probability of the observed choice
  if (!is.null(use)) ll <- ll[use]
  -sum(ll)
}

# Log-likelihood of a piece's observed choices at every sigma of the grid, as a
# vector over the grid (only the trials where `use` is TRUE; NULL = all). One
# boxes x grid computation per piece.
piece_loglik_grid <- function(piece, variant, par, use = NULL, exact = TRUE) {
  sigma <- matrix(SIGMA_GRID, nrow(piece$boxes), length(SIGMA_GRID), byrow = TRUE)
  lp    <- piece_choice_logprob(piece, variant, sigma, par, exact)   # trials x grid, log probabilities of left and of right
  ll    <- lp$right
  isL   <- piece$trials$ChooseL == 1
  ll[isL, ] <- lp$left[isL, ]                                  # the observed choice's log probability
  if (!is.null(use)) ll <- ll[use, , drop = FALSE]
  colSums(ll)
}

# The train mask of a piece, in that piece's trial order (NULL = all trials).
mask_for <- function(piece, train) if (is.null(train)) NULL else train[piece$trials$global_id]

# Subject step: the log-likelihood of every subject on the sigma grid given the
# group parameters, a subjects x grid matrix with rows named by Subject. The
# subjects are independent given those parameters, so one pass is exact. A
# subject with no training trials gets a row of zeros (no information: their
# posterior is the population distribution). exact = FALSE uses the probit for
# every trial (the approximation the exact calculation replaced; checks only).
grid_pass <- function(prep, par, train = NULL, subjects = prep$subjects$Subject, exact = TRUE) {
  ll <- matrix(0, length(subjects), length(SIGMA_GRID), dimnames = list(subjects, NULL))
  for (s in subjects) {
    piece <- prep$subj[[s]]
    use   <- mask_for(piece, train)
    if (!is.null(use) && !any(use)) next
    ll[s, ] <- piece_loglik_grid(piece, prep$variant, par, use, exact)
  }
  ll
}

# Log weights of the population distribution on the grid: a normal on log
# sigma with mean mu and SD tau, truncated to the grid (the weights sum to 1).
population_logweights <- function(mu, tau) {
  lw <- dnorm(SIGMA_GRID_LOG, mu, tau, log = TRUE)
  m  <- max(lw)
  lw - m - log(sum(exp(lw - m)))
}

# Marginal negative log-likelihood of a set of subjects given (mu, tau):
# -sum_i log sum_k exp(ll_ik + logw_k). Takes the rows' maxima rmax and
# E = exp(ll - rmax), computed once per grid pass, so that nothing overflows
# and each evaluation is one matrix-vector product.
marginal_nll <- function(E, rmax, mu, tau) {
  w <- exp(population_logweights(mu, tau))
  -sum(rmax + log(as.vector(E %*% w)))
}

# Population step: the (mu, tau) of one range condition that maximise the
# marginal likelihood of its subjects, given their grid log-likelihoods `ll`
# (subjects x grid). Nelder-Mead on the logistic-bounded scale, started from
# the median and spread of the subjects' own best grid points. Returns
# list(mu, tau, nll, convergence, post), where post (subjects x grid) is each
# subject's posterior over sigma given their trials and the fitted population.
fit_population <- function(ll) {
  rmax <- apply(ll, 1, max)
  E    <- exp(ll - rmax)
  best <- SIGMA_GRID_LOG[apply(ll, 1, which.max)]
  mu0  <- clamp(median(best), MU_BOUNDS[1], MU_BOUNDS[2])
  tau0 <- clamp(IQR(best) / 1.349, TAU_BOUNDS[1], TAU_BOUNDS[2])
  f    <- function(p) marginal_nll(E, rmax, mu_th(p[1]), tau_th(p[2]))
  opt  <- optim(c(th_mu(mu0), th_tau(tau0)), f, method = "Nelder-Mead", control = INNER_CONTROL)
  mu   <- mu_th(opt$par[1])
  tau  <- tau_th(opt$par[2])
  post <- sweep(E, 2, exp(population_logweights(mu, tau)), "*")
  post <- post / rowSums(post)
  list(mu = mu, tau = tau, nll = opt$value, convergence = opt$convergence, post = post)
}

# The worst convergence code of the population steps of the current fit
# (0 = every Nelder-Mead stopped on its tolerance). Reset by fit_profile.
.inner <- new.env()
.inner$worst <- 0L

# Helpers to move one condition's parameters in and out of `par`.
set_cond <- function(par, cc, s0 = NULL, q = NULL, lambda = NULL) {
  if (!is.null(s0))     par$s0[cc]     <- s0
  if (!is.null(q))      par$q[cc]      <- q
  if (!is.null(lambda)) par$lambda[cc] <- lambda
  par
}
new_par <- function(conds, s0 = 1, q = 0.5, lambda = 1) {
  list(s0 = setNames(rep(s0, length(conds)), conds), q = setNames(rep(q, length(conds)), conds),
       lambda = setNames(rep(lambda, length(conds)), conds))
}

# Which group parameters a variant fits besides the population of sigma. q is
# free only for a fitted prior: the zero and flat priors fix it at 0.5, the
# middle of the span, which is an objective 0 for the symmetric spans of this
# design. s0 is 1, or 1e3 for the flat prior.
free_q_of   <- function(variant) variant$prior == "fitted"
s0_fixed_of <- function(variant) if (variant$prior == "flat") 1e3 else 1

# Profile NLL of one range condition at `par`: its subjects' grid
# log-likelihoods (or the rows of an existing grid pass `ll`), then the
# population step. Returns the list of fit_population.
profile_cond <- function(prep, cc, par, train = NULL, ll = NULL) {
  subs <- prep$subjects$Subject[prep$subjects$WN == cc]
  if (is.null(ll)) ll <- grid_pass(prep, par, train, subs)
  pop <- fit_population(ll[subs, , drop = FALSE])
  .inner$worst <- max(.inner$worst, pop$convergence)
  pop
}

# Profile NLL over both conditions at `par`: one grid pass, the population
# step per condition, the NLLs summed.
profile_total <- function(prep, par, train = NULL, exact = TRUE) {
  ll <- grid_pass(prep, par, train, exact = exact)
  total <- 0
  for (cc in names(prep$cond)) total <- total + profile_cond(prep, cc, par, train, ll)$nll
  total
}

# Fit on the profile likelihood, in the three stages of the header, to all
# trials or (with a `train` mask) to the training trials only. Returns
# list(sigma, mu, tau, s0, q, lambda, nll, grid, variant, post, convergence,
# convergence_inner): sigma is each subject's posterior mean, post their
# posterior over the sigma grid (subjects x grid), convergence holds the code
# of every outer Nelder-Mead call (0 = converged within joint_maxit function
# evaluations; empty when none was needed) and convergence_inner the worst
# code of the population steps.
fit_profile <- function(prep, q_grid = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                        joint_maxit = 200, train = NULL, verbose = TRUE) {
  variant  <- prep$variant
  conds    <- names(prep$cond)
  free_q   <- free_q_of(variant)
  lam_mode <- variant$lambda
  s0_fixed <- s0_fixed_of(variant)
  conv     <- c()
  .inner$worst <- 0L
  if (!free_q) q_grid <- 0.5

  # Stage 1: grid over q with lambda = 1. One grid pass per point serves both
  # conditions (a single point, q = 0.5, when q is fixed).
  grid <- tibble()
  for (q in q_grid) {
    par <- new_par(conds, s0_fixed, q, 1)
    ll  <- grid_pass(prep, par, train)
    row <- c()
    for (cc in conds) {
      pop     <- profile_cond(prep, cc, par, train, ll)
      row[cc] <- pop$nll
      grid    <- bind_rows(grid, tibble(cond = cc, q = q, nll = pop$nll, mu = pop$mu, tau = pop$tau))
    }
    if (verbose) cat(sprintf("  grid q = %5.2f  nll = %s\n", q, paste(sprintf("%.1f", row), collapse = " / ")))
  }

  # Stage 2: refine each condition from its best grid point. Lambda is only
  # refined here when it is per condition; otherwise it stays at 1 for now.
  # Bounded Nelder-Mead for two parameters, a bounded 1-D search for one.
  par <- new_par(conds, s0_fixed)
  for (cc in conds) {
    best <- grid %>% filter(cond == cc) %>% slice_min(nll, n = 1, with_ties = FALSE)
    b_q  <- clamp(best$q, Q_BOUNDS[1], Q_BOUNDS[2])
    if (!free_q && lam_mode == "cond") {
      f   <- function(th) profile_cond(prep, cc, set_cond(par, cc, q = 0.5, lambda = lam_th(th)), train)$nll
      opt <- optimize(f, interval = th_lam(c(LAMBDA_BOUNDS[1] * 1.01, LAMBDA_BOUNDS[2] * 0.99)))
      par <- set_cond(par, cc, q = 0.5, lambda = lam_th(opt$minimum))
    } else if (!free_q) {
      par <- set_cond(par, cc, q = 0.5, lambda = 1)
    } else if (lam_mode == "cond") {
      f   <- function(p) profile_cond(prep, cc, set_cond(par, cc, q = q_th(p[1]), lambda = lam_th(p[2])), train)$nll
      opt <- optim(c(th_q(b_q), th_lam(1)), f, method = "Nelder-Mead", control = list(maxit = 100, reltol = 1e-7))
      par <- set_cond(par, cc, q = q_th(opt$par[1]), lambda = lam_th(opt$par[2]))
      conv[cc] <- opt$convergence
    } else {
      f   <- function(q) profile_cond(prep, cc, set_cond(par, cc, q = q, lambda = 1), train)$nll
      opt <- optimize(f, interval = c(max(Q_BOUNDS[1], b_q - 0.3), min(Q_BOUNDS[2], b_q + 0.3)))
      par <- set_cond(par, cc, q = opt$minimum, lambda = 1)
    }
    if (verbose) cat(sprintf("  %s: q = %.3f  lambda = %.3f\n", cc, par$q[[cc]], par$lambda[[cc]]))
  }

  # Stage 3: one lambda shared by both conditions, refined jointly with the
  # two q's from the per-condition solutions on the bounded scale; with q
  # fixed, a bounded 1-D search over lambda alone.
  if (lam_mode == "shared") {
    if (free_q) {
      unpack <- function(p) {
        out <- par
        for (i in seq_along(conds)) out$q[conds[i]] <- q_th(p[i])
        out$lambda[] <- lam_th(p[length(p)])
        out
      }
      f   <- function(p) profile_total(prep, unpack(p), train)
      p0  <- c(th_q(clamp(unname(par$q[conds]), Q_BOUNDS[1], Q_BOUNDS[2])), th_lam(1))
      # For Nelder-Mead, optim's maxit caps the number of function evaluations.
      opt <- optim(p0, f, method = "Nelder-Mead", control = list(maxit = joint_maxit, reltol = 1e-7))
      par <- unpack(opt$par)
      conv["joint"] <- opt$convergence
    } else {
      f   <- function(th) { pp <- par; pp$lambda[] <- lam_th(th); profile_total(prep, pp, train) }
      opt <- optimize(f, interval = th_lam(c(LAMBDA_BOUNDS[1] * 1.01, LAMBDA_BOUNDS[2] * 0.99)))
      par$lambda[] <- lam_th(opt$minimum)
    }
    if (verbose) cat(sprintf("  shared lambda = %.3f  q = %s\n", par$lambda[[1]], paste(sprintf("%.3f", par$q), collapse = " / ")))
  }

  # The population and every subject's posterior at the fitted group parameters.
  ll   <- grid_pass(prep, par, train)
  nll  <- 0
  mu   <- c()
  tau  <- c()
  post <- ll * 0
  for (cc in conds) {
    pop     <- profile_cond(prep, cc, par, train, ll)
    nll     <- nll + pop$nll
    mu[cc]  <- pop$mu
    tau[cc] <- pop$tau
    post[rownames(pop$post), ] <- pop$post
  }
  post  <- post[prep$subjects$Subject, , drop = FALSE]
  sigma <- setNames(as.vector(post %*% SIGMA_GRID), rownames(post))
  if (verbose) cat(sprintf("  population: exp(mu) = %s  tau = %s  nll = %.2f\n",
                           paste(sprintf("%.2f", exp(mu)), collapse = " / "), paste(sprintf("%.2f", tau), collapse = " / "), nll))
  list(sigma = sigma, mu = mu, tau = tau, s0 = par$s0, q = par$q, lambda = par$lambda, nll = nll, grid = grid,
       variant = variant, post = post, convergence = conv, convergence_inner = .inner$worst)
}

# P(choose left) for every trial, averaged over each subject's posterior on
# the sigma grid (`post`, subjects x grid, from fit_profile): the model's
# prediction for a trial given what the subject's fitted trials say about
# their noise. Returns prep$all$trials plus p_left. `trials` (a logical over
# prep$all$trials) restricts the work to those trials; the others get NA.
posterior_choice_probabilities <- function(prep, post, par, trials = NULL) {
  out <- prep$all$trials
  out$p_left <- NA_real_
  for (s in rownames(post)) {
    piece <- prep$subj[[s]]
    keep  <- mask_for(piece, trials)
    if (!is.null(keep) && !any(keep)) next
    sigma <- matrix(SIGMA_GRID, nrow(piece$boxes), length(SIGMA_GRID), byrow = TRUE)
    p     <- as.vector(piece_choice_prob(piece, prep$variant, sigma, par) %*% post[s, ])
    if (!is.null(keep)) p[!keep] <- NA_real_
    out$p_left[piece$trials$global_id] <- p
  }
  out
}

# Cross-validation over trials. Each subject's trials are numbered in Trial
# order and cut into blocks of TRIALS_PER_FOLD; fold k holds out block k of
# every subject. A subject whose trial count is not a multiple of
# TRIALS_PER_FOLD has a shorter last block, and a subject with fewer than
# (k - 1) * TRIALS_PER_FOLD + 1 trials has nothing held out in fold k and is
# entirely in training (here folds 7-12 hold out trials of the 60-trial
# subjects only). Every trial is held out exactly once. The model is refit
# from scratch on everything else by
# fit_profile (the same three stages as the full fit, with the population
# distribution refit on the training trials) and predicts the held-out trials
# from each subject's posterior over sigma given their training trials. Returns
#   $pred         prep$all$trials plus p_left_heldout and fold
#   $fold_params  s0, q, lambda, mu and tau fit in every fold, the training
#                 NLL, and the worst convergence codes of the fold (0 = all
#                 converged), outer and inner
# Cached to cache_file (and returned from it when it exists).
TRIALS_PER_FOLD <- 5
loo_by_trial_index <- function(prep, cache_file, verbose = TRUE, folds = NULL) {
  if (file.exists(cache_file)) return(readRDS(cache_file))
  trials <- prep$all$trials %>%
    group_by(Subject) %>% mutate(k = rank(Trial)) %>% ungroup() %>%
    arrange(global_id)
  K <- ceiling(max(trials$k) / TRIALS_PER_FOLD)
  # Every fold is saved to its own file, so the folds can be split over
  # several processes (`folds` = the ones this process does; NULL = all) and
  # an interrupted run resumes. The combined cache is written once all K exist.
  fold_dir  <- sub("\\.rds$", "_folds", cache_file)
  fold_file <- function(k) file.path(fold_dir, sprintf("fold_%02d.rds", k))
  dir.create(fold_dir, showWarnings = FALSE, recursive = TRUE)
  if (is.null(folds)) folds <- seq_len(K)
  for (k in folds) {
    if (file.exists(fold_file(k))) next
    held  <- ceiling(trials$k / TRIALS_PER_FOLD) == k
    train <- !held
    t0  <- Sys.time()
    fit <- fit_profile(prep, train = train, verbose = FALSE)
    pr  <- posterior_choice_probabilities(prep, fit$post, list(s0 = fit$s0, q = fit$q, lambda = fit$lambda), trials = held)
    stopifnot(!anyNA(pr$p_left[held]))
    conv <- if (length(fit$convergence) > 0) max(fit$convergence) else 0L
    params <- data.frame(fold = k, cond = names(fit$s0), s0 = unname(fit$s0),
                         q = unname(fit$q[names(fit$s0)]), lambda = unname(fit$lambda[names(fit$s0)]),
                         mu = unname(fit$mu[names(fit$s0)]), tau = unname(fit$tau[names(fit$s0)]),
                         nll_train = fit$nll, convergence = conv, convergence_inner = fit$convergence_inner)
    saveRDS(list(fold = k, global_id = trials$global_id[held], p_left = pr$p_left[held], params = params), fold_file(k))
    if (verbose) cat(sprintf("fold %d of %d: %d held-out trials  q = %s  lambda = %s  exp(mu) = %s  tau = %s  train nll = %.1f  convergence = %d / %d  %.1f min  %s\n",
                             k, K, sum(held), paste(sprintf("%.3f", fit$q), collapse = " / "),
                             paste(sprintf("%.3f", fit$lambda), collapse = " / "),
                             paste(sprintf("%.2f", exp(fit$mu)), collapse = " / "), paste(sprintf("%.2f", fit$tau), collapse = " / "),
                             fit$nll, conv, fit$convergence_inner,
                             as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time(), "%H:%M:%S")))
  }
  done <- file.exists(sapply(seq_len(K), fold_file))
  if (!all(done)) {
    if (verbose) cat(sprintf("%d of %d folds done; the combined file is written when all exist\n", sum(done), K))
    return(NULL)
  }
  p_heldout <- rep(NA_real_, nrow(trials))
  fold_of   <- rep(NA_integer_, nrow(trials))
  fold_params <- list()
  for (k in seq_len(K)) {
    f   <- readRDS(fold_file(k))
    idx <- match(f$global_id, trials$global_id)
    p_heldout[idx]   <- f$p_left
    fold_of[idx]     <- k
    fold_params[[k]] <- f$params
  }
  out <- list(
    pred        = trials %>% mutate(p_left_heldout = p_heldout, fold = fold_of) %>% select(-k),
    fold_params = bind_rows(fold_params),
    variant     = prep$variant
  )
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(out, cache_file)
  out
}

# Out-of-sample log-likelihood of held-out predictions, per trial.
heldout_loglik <- function(pred) {
  p <- pmin(pmax(pred$p_left_heldout, 1e-12), 1 - 1e-12)
  ifelse(pred$ChooseL == 1, log(p), log(1 - p))
}

# ---------------------------------------------------------------------------
# The model variants, and one function that fits a variant and runs its
# leave-one-out with caching (used by code/fitting/run_bayesian_fits.R from a
# terminal and by notebooks/bayesian_model.Rmd, which then only loads caches).
# ---------------------------------------------------------------------------

# units / prior_ref: which reference distribution sets the units of value and
# noise, and which the prior is calibrated to ("type" = per attribute type x
# range condition; "cond" = per range condition; "points" = no rescaling).
# Because every subject is in one range condition and has their own noise, the
# units only matter within a subject, for how the attribute types are weighted
# against each other; the current variants all use points (equal weighting).
# prior: "fitted" (position q per condition), "zero" (an objective 0, q fixed
# at 0.5) or "flat" (mean at 0, SD 1e3 x the reference spread: no pull toward
# the prior, but the loss-aversion kink still sits at 0).
# lambda: "shared" (one loss-aversion coefficient for both conditions), "cond"
# (one per condition) or "fixed1" (no loss aversion). The same algorithm
# (fit_profile) fits every variant.
bayesian_variants <- function() {
  main <- function(prior_ref, prior, lambda, label)
    list(units = "points", prior_ref = prior_ref, prior = prior, lambda = lambda, label = label)
  list(
    cond_lam     = main("cond", "fitted", "shared", "Per-condition prior, loss aversion"),
    cond_lam1    = main("cond", "fitted", "fixed1", "Per-condition prior"),
    type_lam     = main("type", "fitted", "shared", "Per-type prior, loss aversion"),
    type_lam1    = main("type", "fitted", "fixed1", "Per-type prior"),
    cond_lamcond = main("cond", "fitted", "cond",   "Per-condition prior, loss aversion per range"),
    zero         = main("cond", "zero",   "shared", "Zero prior, loss aversion"),
    flat         = main("cond", "flat",   "shared", "Flat prior, loss aversion")
  )
}

# Full-data fit (fit_<name>.rds) and leave-one-out (loo_<name>.rds) of one
# variant, both read from cache_dir when they exist. Returns
# list(name, variant, prep, fit, loo).
run_variant <- function(name, boxlong, cache_dir, verbose = TRUE, what = c("both", "fit", "loo"), folds = NULL) {
  what <- match.arg(what)
  v    <- bayesian_variants()[[name]]
  prep <- prepare_boxes(boxlong, v)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  fit_file <- file.path(cache_dir, paste0("fit_", name, ".rds"))
  fit <- NULL
  if (what != "loo") {
    if (file.exists(fit_file)) {
      fit <- readRDS(fit_file)
    } else {
      if (verbose) cat("\n==== fitting", name, " ", format(Sys.time()), "\n")
      t0  <- Sys.time()
      fit <- fit_profile(prep, verbose = verbose)
      if (verbose) cat(sprintf("  convergence: %s  inner: %d  %.1f min\n", paste(fit$convergence, collapse = " "),
                               fit$convergence_inner, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      saveRDS(fit, fit_file)
    }
  }
  loo <- NULL
  if (what != "fit") {
    if (verbose) cat("\n==== cross-validation", name, " ", format(Sys.time()), "\n")
    loo <- loo_by_trial_index(prep, file.path(cache_dir, paste0("loo_", name, ".rds")), verbose = verbose, folds = folds)
  }
  list(name = name, variant = v, prep = prep, fit = fit, loo = loo)
}
