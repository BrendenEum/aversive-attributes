# # BRMS default settings
# brm <- function(...)
#   brms::brm(
#     ...,
#     iter = 6000,
#     warmup = floor(6000/2),
#     chains = 4,
#     cores = 4,
#     seed = 4,
#     backend = "cmdstanr",
#     silent = 0,
#     refresh = 6000/50,
#     file_refit = "on_change"
#   )
