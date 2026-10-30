# Maximum-likelihood fitting of the Bayesian attribute-level model (the model
# itself is in code/simulation/simulate_bayesian_attributes.R; source that
# first). Everything uses the closed-form choice probabilities, so no
# simulation is needed to fit.
#
# Parameters
#   sigma   noise per second, one per subject (pinned down by the subject's accuracy)
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
# Error model. The sigmas are partially pooled (variant$shrink = TRUE): each
# subject's log sigma gets a Gaussian pull toward the group median with SD
# SIGMA_SHRINK_TAU, the same logic as a subject random effect. It only bites
# for subjects whose own trials leave sigma poorly determined (e.g. 30 correct
# trials, for which plain maximum likelihood says sigma = 0 and a held-out
# error then costs 25+ log-likelihood units). The reported NLL is the plain
# likelihood, without the penalty.
#
# One algorithm fits every variant; the variants differ only in which group
# parameters are free. Given the group parameters, every subject's sigma is a
# separate 1-D problem, so the fit works on the profile likelihood
# (fit_profile), in three stages: (1) for a grid of q with lambda = 1 the
# sigmas are re-optimised and the NLL recorded; (2) q is refined per condition
# from the best grid point by a bounded 1-D search (together with a
# per-condition lambda, by Nelder-Mead, when lambda = "cond"); (3) for a shared
# lambda, q in both conditions and lambda are refined jointly by Nelder-Mead.
# The sigmas are re-optimised at every candidate. When q is fixed (zero and
# flat priors) stage 1 is a single sigma pass, stage 2 has nothing to refine
# and stage 3 is a bounded 1-D search over lambda. The leave-one-out
# (loo_by_trial_index) fits every fold from scratch by the same three stages
# on the training trials, shrinkage centre included; nothing is carried over
# from the full fit.
#
# Variants: prior = "zero" fixes the prior position at an objective 0 (q = 0.5,
# the middle of the symmetric span); prior = "flat" keeps that position and
# sets s0 = 1e3, so the prior carries no weight in the estimate but the kink of
# the loss aversion still sits at 0. Both fit the sigmas and lambda.
#
# Group parameters travel in a list `par` with s0, q, lambda, each a vector
# named by condition. Train/test masks are logical vectors over
# prep$all$trials (global_id order).

SIGMA_BOUNDS  <- c(0.005, 500)   # search range for sigma (log scale), wide enough for any units
Q_BOUNDS      <- c(-1, 2)        # q may sit outside the reference range
LAMBDA_BOUNDS <- c(0.25, 8)      # loss aversion (log scale); below 1 = gain seeking
SIGMA_SHRINK_TAU <- 1            # SD (log units) of the pull of each log sigma toward the group median

# Nelder-Mead has no bounds, so the group parameters are optimised on an
# unconstrained scale and mapped into their bounds with a logistic: without
# this, a flat likelihood (e.g. a weak prior, where q no longer matters) lets
# the simplex drift to absurd values.
to_theta   <- function(x, lo, hi) { p <- (x - lo) / (hi - lo); log(p / (1 - p)) }
from_theta <- function(th, lo, hi) lo + (hi - lo) / (1 + exp(-th))
th_q   <- function(q)   to_theta(q, Q_BOUNDS[1], Q_BOUNDS[2])
q_th   <- function(th)  from_theta(th, Q_BOUNDS[1], Q_BOUNDS[2])
th_lam <- function(lam) to_theta(log(lam), log(LAMBDA_BOUNDS[1]), log(LAMBDA_BOUNDS[2]))
lam_th <- function(th)  exp(from_theta(th, log(LAMBDA_BOUNDS[1]), log(LAMBDA_BOUNDS[2])))
clamp  <- function(x, lo, hi) pmin(pmax(x, lo + 1e-6 * (hi - lo)), hi - 1e-6 * (hi - lo))

# Negative log-likelihood of a piece's observed choices (only the trials where
# `use` is TRUE; NULL = all). sigma is a per-box vector or one number.
piece_nll <- function(piece, variant, sigma, par, use = NULL) {
  mm <- piece_moments(piece, variant, sigma, par)
  zs <- ifelse(mm$var > 0, mm$delta / sqrt(mm$var), 0)
  # Log of the normal CDF taken directly so confident trials don't underflow.
  ll <- ifelse(piece$trials$ChooseL == 1, pnorm(zs, log.p = TRUE), pnorm(-zs, log.p = TRUE))
  if (!is.null(use)) ll <- ll[use]
  -sum(ll)
}

# The train mask of a piece, in that piece's trial order (NULL = all trials).
mask_for <- function(piece, train) if (is.null(train)) NULL else train[piece$trials$global_id]

# Subject step: the sigma of one subject, given the group parameters. With
# shrinkage, the objective adds (log sigma - center)^2 / (2 tau^2), where
# `center` is the group median of log sigma from an unpenalised pass.
fit_subject_sigma <- function(piece, variant, par, use = NULL, center = NULL) {
  shrink <- isTRUE(variant$shrink) && !is.null(center)
  f <- function(log_sigma) {
    piece_nll(piece, variant, exp(log_sigma), par, use) +
      (if (shrink) (log_sigma - center)^2 / (2 * SIGMA_SHRINK_TAU^2) else 0)
  }
  opt <- optimize(f, interval = log(SIGMA_BOUNDS))
  exp(opt$minimum)
}

# The sigmas of a set of subjects given the group parameters: one pass is the
# exact conditional optimum because the sigmas are independent given those
# parameters. Returns a vector named by Subject (NA for a subject with no
# training trials).
sigma_pass <- function(prep, par, train = NULL, subjects = prep$subjects$Subject) {
  sigma <- setNames(rep(NA_real_, length(subjects)), subjects)
  for (s in subjects) {
    piece <- prep$subj[[s]]
    use   <- mask_for(piece, train)
    if (!is.null(use) && !any(use)) next
    sigma[s] <- fit_subject_sigma(piece, prep$variant, par, use, prep$sigma_center)
  }
  sigma
}

# The centre of the shrinkage: the group median of log sigma from one
# unpenalised pass at the starting parameters (over the training trials when
# `train` is given). Stored in prep$sigma_center and carried in the fit object.
sigma_shrink_center <- function(prep, train = NULL) {
  if (!isTRUE(prep$variant$shrink)) return(NULL)
  prep$sigma_center <- NULL
  sigma <- sigma_pass(prep, new_par(names(prep$cond)), train)
  median(log(sigma), na.rm = TRUE)
}

# NLL of one condition given the sigmas (a vector named by Subject).
cond_nll <- function(prep, cc, sigma, par, train = NULL) {
  piece <- prep$cond[[cc]]
  piece_nll(piece, prep$variant, unname(sigma[as.character(piece$boxes$Subject)]), par, mask_for(piece, train))
}

# Profile NLL of one range condition at `par`: that condition's sigmas
# re-optimised, then its NLL. Returns list(nll, sigma).
profile_cond <- function(prep, cc, par, train = NULL) {
  subs  <- prep$subjects$Subject[prep$subjects$WN == cc]
  sigma <- sigma_pass(prep, par, train, subs)
  list(nll = cond_nll(prep, cc, sigma, par, train), sigma = sigma)
}

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

# Which group parameters a variant fits besides the sigmas. q is free only for
# a fitted prior: the zero and flat priors fix it at 0.5, the middle of the
# span, which is an objective 0 for the symmetric spans of this design. s0 is
# 1, or 1e3 for the flat prior.
free_q_of   <- function(variant) variant$prior == "fitted"
s0_fixed_of <- function(variant) if (variant$prior == "flat") 1e3 else 1

# Profile NLL over both conditions at `par`: every sigma re-optimised, then
# the NLLs summed.
profile_total <- function(prep, par, train = NULL) {
  sigma <- sigma_pass(prep, par, train)
  total <- 0
  for (cc in names(prep$cond)) total <- total + cond_nll(prep, cc, sigma, par, train)
  total
}

# Fit on the profile likelihood, in the three stages of the header, to all
# trials or (with a `train` mask) to the training trials only. Returns
# list(sigma, s0, q, lambda, nll, grid, variant, sigma_center, convergence),
# where convergence holds the code of every Nelder-Mead call (0 = converged
# within joint_maxit function evaluations); empty when none was needed.
fit_profile <- function(prep, q_grid = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7),
                        joint_maxit = 200, train = NULL, verbose = TRUE) {
  variant  <- prep$variant
  conds    <- names(prep$cond)
  free_q   <- free_q_of(variant)
  lam_mode <- variant$lambda
  s0_fixed <- s0_fixed_of(variant)
  conv     <- c()
  if (!free_q) q_grid <- 0.5
  prep$sigma_center <- sigma_shrink_center(prep, train)
  if (verbose && !is.null(prep$sigma_center)) cat(sprintf("  shrinkage centre: median sigma = %.3f\n", exp(prep$sigma_center)))

  # Stage 1: grid over q with lambda = 1. One sigma pass per point serves both
  # conditions (a single point, q = 0.5, when q is fixed).
  grid <- tibble()
  for (q in q_grid) {
    par   <- new_par(conds, s0_fixed, q, 1)
    sigma <- sigma_pass(prep, par, train)
    row   <- c()
    for (cc in conds) {
      row[cc] <- cond_nll(prep, cc, sigma, par, train)
      grid    <- bind_rows(grid, tibble(cond = cc, q = q, nll = row[cc]))
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

  sigma <- sigma_pass(prep, par, train)
  nll   <- 0
  for (cc in conds) nll <- nll + cond_nll(prep, cc, sigma, par, train)
  if (verbose) cat(sprintf("  nll = %.2f\n", nll))
  list(sigma = sigma[prep$subjects$Subject], s0 = par$s0, q = par$q, lambda = par$lambda, nll = nll, grid = grid,
       variant = variant, sigma_center = prep$sigma_center, convergence = conv)
}

# Leave-one-out over trials. Fold k holds out the k-th trial (in Trial order)
# of every subject; the model is refit from scratch on everything else by
# fit_profile (the same three stages as the full fit, with the shrinkage
# centre recomputed from the training trials) and predicts the held-out
# trials. Returns
#   $pred         prep$all$trials plus p_left_heldout and fold
#   $fold_params  s0, q and lambda fit in every fold, the training NLL, and
#                 the worst Nelder-Mead convergence code of the fold (0 = all converged)
# Cached to cache_file (and returned from it when it exists).
loo_by_trial_index <- function(prep, cache_file, verbose = TRUE) {
  if (file.exists(cache_file)) return(readRDS(cache_file))
  trials <- prep$all$trials %>%
    group_by(Subject) %>% mutate(k = rank(Trial)) %>% ungroup() %>%
    arrange(global_id)
  K <- max(trials$k)
  p_heldout <- rep(NA_real_, nrow(trials))
  fold_of   <- rep(NA_integer_, nrow(trials))
  fold_params <- list()
  for (k in seq_len(K)) {
    held  <- trials$k == k
    train <- !held
    t0  <- Sys.time()
    fit <- fit_profile(prep, train = train, verbose = FALSE)
    stopifnot(!anyNA(fit$sigma))   # every subject keeps all but one trial, so every sigma is fit
    pr  <- choice_probabilities(prep, fit$sigma, list(s0 = fit$s0, q = fit$q, lambda = fit$lambda))
    p_heldout[held] <- pr$p_left[held]
    fold_of[held]   <- k
    conv <- if (length(fit$convergence) > 0) max(fit$convergence) else 0L
    fold_params[[k]] <- data.frame(fold = k, cond = names(fit$s0), s0 = unname(fit$s0),
                                   q = unname(fit$q[names(fit$s0)]), lambda = unname(fit$lambda[names(fit$s0)]),
                                   nll_train = fit$nll, convergence = conv)
    if (verbose) cat(sprintf("fold %d of %d: %d held-out trials  q = %s  lambda = %s  train nll = %.1f  convergence = %d  %.1f min  %s\n",
                             k, K, sum(held), paste(sprintf("%.3f", fit$q), collapse = " / "),
                             paste(sprintf("%.3f", fit$lambda), collapse = " / "), fit$nll, conv,
                             as.numeric(difftime(Sys.time(), t0, units = "mins")), format(Sys.time(), "%H:%M:%S")))
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
# (one per condition) or "fixed1" (no loss aversion). shrink: partial pooling
# of the sigmas. The same algorithm (fit_profile) fits every variant.
bayesian_variants <- function() {
  main <- function(prior_ref, prior, lambda, label)
    list(units = "points", prior_ref = prior_ref, prior = prior, lambda = lambda, shrink = TRUE, label = label)
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
run_variant <- function(name, boxlong, cache_dir, verbose = TRUE) {
  v    <- bayesian_variants()[[name]]
  prep <- prepare_boxes(boxlong, v)
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  fit_file <- file.path(cache_dir, paste0("fit_", name, ".rds"))
  if (file.exists(fit_file)) {
    fit <- readRDS(fit_file)
  } else {
    if (verbose) cat("\n==== fitting", name, " ", format(Sys.time()), "\n")
    t0  <- Sys.time()
    fit <- fit_profile(prep, verbose = verbose)
    if (verbose) cat(sprintf("  convergence: %s  %.1f min\n", paste(fit$convergence, collapse = " "),
                             as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    saveRDS(fit, fit_file)
  }
  if (verbose) cat("\n==== leave-one-out", name, " ", format(Sys.time()), "\n")
  loo <- loo_by_trial_index(prep, file.path(cache_dir, paste0("loo_", name, ".rds")), verbose = verbose)
  list(name = name, variant = v, prep = prep, fit = fit, loo = loo)
}
