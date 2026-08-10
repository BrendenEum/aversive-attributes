# function to simulate DDM

simulate_ddm <- function(evidence, a = 1.0, d = 0.01, s = 0.02, b = 0.0, c = 0.0, ndt = 1, ts_size = .01, rt_max = 60) {
  
  rt_max = rt_max*1000
  ts_size = ts_size*1000
  
  bounds = a
  rdv_t = b
  inside_bounds = T
  timestep = 0
  
  while(inside_bounds==T & timestep<(rt_max/ts_size)) {
    rdv_t = rdv_t + d * evidence + rnorm(1, 0, s)
    timestep = timestep + 1
    if (abs(rdv_t)>=bounds) {inside_bounds = F}
  }
  
  if (rdv_t >= bounds) {choice = 1}
  if (rdv_t <= bounds) {choice = 0}
  rt = (timestep * ts_size) + (ndt*1000)
  if (rt >= rt_max) {choice = NA}
  
  return(list(choice=choice, rt=rt/1000))
  
}