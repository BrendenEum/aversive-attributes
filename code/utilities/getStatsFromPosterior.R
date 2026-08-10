getStatsFromPosterior = function(posterior_samples) {
  b = mean(posterior_samples) %>% round(3)
  b_se = sd(posterior_samples) %>% round(3)
  b_ci = quantile(posterior_samples, c(.025, .975))
  output = list(mean = b, se = b_se, conf_int = b_ci)
  return(output)
}