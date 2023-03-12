#' Title
#'
#' @param n_trts 
#' @param n_blocks 
#' @param period_length 
#' @param washout_length 
#' @param sampling_timestep 
#' @param sd_b 
#' @param sd_p 
#' @param sd_o 
#' @param Ejs 
#' @param gammas 
#' @param taus 
#' @param eta 
#' @param custom.Bt 
#' @param Xj_inits 
#' @param Z_init 
#'
#' @return
#' @export simulate
#'
#' @examples
simulate = function(n_trts, n_blocks, period_length, washout_length,
                    sampling_timestep, sd_b, sd_p, sd_o,
                    Ejs, gammas, taus, eta, 
                    custom.Bt = NULL, Xj_inits = NULL, Z_init = NULL) {
  
  # TO-DO:
  # incorporate an interaction effects?
  # 
  
  # Establish timescales over simulated trial
  duration = n_blocks * n_trts * (period_length + washout_length) - washout_length
  timesteps = seq(0, duration, sampling_timestep)
  time_init = timesteps[1]
  time_final = timesteps[length(timesteps)]
  
  # Adjust number of periods + treatment order based on whether or not washout is used
  if (washout_length == 0) {
    n_periods = n_blocks * n_trts 
    order = combinat::permn(1:n_trts) %>% 
      sample(size = n_blocks, replace = T) %>% unlist()
    switches = seq(0, duration, period_length)
  } else {
    n_periods = (n_blocks * (2 * n_trts)) - 1
    order = combinat::permn(1:n_trts) %>% 
      sample(size = n_blocks, replace = T) %>% unlist() %>% 
      incorporate_washout()
    switches = find_treatment_switches(n_blocks * n_trts, period_length, washout_length)
  }
  
  # Construct necessary functions for simulations
  Bt = simulate_baseline_function(timesteps, sd_b, sampling_timestep)
  
  # Treatment number by timestep
  Tt = approxfun(switches, c(order, order[length(order)]), method = 'constant')
  
  # Construct the effect functions for each treatment
  Xjt = lapply(1:n_trts, function(j) { 
    Tjt = as.numeric(Tt(timesteps) == j)
    simulate_effect_function(j, order, Tjt, switches, Xj_init = 0, Ejs[j], taus[j], gammas[j] )
  })
  
  # Get the total treatment effect across time
  sum_Xj = lapply(Xjt, function(f) { f(timesteps) }) %>% Reduce("+", .)
  
  Zt = simulate_outcome(timesteps, Bt, sum_Xj, Z_init = 0, eta, sampling_timestep)
  Yt = Zt(timesteps) + rnorm(length(timesteps), 0, sd_o)
  
  data = tibble(
    time = timesteps,
    period = match_period_to_timesteps(order, switches),
    treatment = Tt(time),
    Bt = Bt(time),
    sum_Xjt = sum_Xj,
    Zt = Zt(time),
    Yt = Yt
  )
  
  for (j in 1:length(Xjt)) {
    data[[paste0("T", j)]] = as.numeric(Tt(timesteps) == j)
    data[[paste0("X", j)]] = Xjt[[j]](timesteps)
  }
  
  out = list()
  out[["data"]] = data
  out[["functions"]] = list(Bt = Bt, 
                            Tt = Tt,
                            Xjt = Xjt,
                            Zt = Zt)
  out[["params"]] = list(
    "trial" = list(n_trts = n_trts,
                   n_blocks = n_blocks,
                   period_length = period_length,
                   washout_length = washout_length),
    "process" = list(timesteps = timesteps,
                     sampling_timestep = sampling_timestep,
                     sd_b = sd_b,
                     sd_p = sd_p,
                     sd_o = sd_o),
    "treatment" = list(Ejs = Ejs,
                       taus = taus,
                       gammas = gammas,
                       eta = eta)
  )
  
  out
  
}

# Incorporate washout into the treatment order
# Convention: washout period will be indicated by the number n_trts+1
incorporate_washout = function(trt_order) {
  order = c()
  for (i in 1:length(trt_order)) {
    if (i != length(trt_order)) {
      order = c(order, c(trt_order[i], n_trts+1))
    } else {
      order = c(order, trt_order[i])              # don't add washout after end of trial
    }
  }
  
  order
}

find_treatment_switches = function(n_periods, period_length, washout_length) {
  switches = c(0)
  val = 0
  for (i in 1:n_periods) {
    val = val + period_length
    switches = c(switches, val)
    val = val + washout_length
    switches = c(switches, val)
  }
  switches = switches[-length(switches)] # removing last addition of washout length
}

simulate_baseline_function = function(time_vec, sd_b, dt) {
  # discrete-time Brownian motion
  wn = rnorm(length(time_vec) - 1, mean = 0, sd = sd_b * sqrt(dt))
  baseline = c(time_vec[1], cumsum(wn))
  approxfun(time_vec, baseline, method = 'linear')
}

simulate_treatment_function = function(n_trts, n_blocks, trt_changes) {
  # counterbalanced + randomized treatment order
  # add n_blocks at end to enable function interpolation
  block_perms = combinat::permn(1:n_trts)
  block_vec = c(sample(block_perms, size = n_blocks, replace = T) %>% unlist(), n_blocks)
  
  list(
    block_vec = block_vec,
    Tt = approxfun(trt_changes, block_vec, method = 'constant')
  )
}

exponential_decay <- function(X_init, target, const, dt) {
  target + (X_init - target) * exp(-dt / const)
}

match_period_to_timesteps = function(order, switches) {
  out = c()
  last = 0
  for (i in 1:(length(order))) {
    len = switches[i+1] - switches[i]
    out = c(out, rep(i, len))
  }
  
  c(out, out[length(out)])
}

create_countup = function(order, switches) {
  # look at t_obs_vec - t_change_vec[period_indices]
  # remember to add the washout days to the previous countup?
  out = c()
  last = 0
  for (i in 1:(length(order))) {
    len = switches[i+1] - switches[i]
    if (i != length(order)) {
      out = c(out, 0:(len-1))
    } else {
      out = c(out, 0:(len))
    }
  }
  
  out
}

simulate_effect_function = function(j, order, Tjt, switches, Xj_init = 0, Ej, tau, gamma) {
  
  starts = switches[1:(length(switches)-1)]
  ends = switches[2:(length(switches))]
  
  # Calculate Xj at the endpoints to figure out if full effect is decayed
  Xj_endpoints = c()
  Xj_endpoints = c(Xj_init)
  n_periods = length(order)
  relevant_periods = which(order == j)
  
  for (period in 1:(n_periods)) {
    
    Tj_used = period %in% relevant_periods
    target = if (Tj_used) Ej else 0 # decide what we're decaying to
    const = if (Tj_used) tau else gamma # decide time constant
    diff = ends[period] - starts[period]
    
    endpoint_val = exponential_decay(Xj_endpoints[period], target, const, diff)
    Xj_endpoints = c(Xj_endpoints, endpoint_val)
    
  }
  
  out = function(timesteps) {
    
    # Washout periods are not new treatments, so they go into the decay
    period_by_timesteps = match_period_to_timesteps(order, switches)
    countup_by_timesteps = create_countup(order, switches)
    target_by_timesteps = ifelse(Tjt, Ej, 0) 
    const_by_timesteps = ifelse(Tjt, tau, gamma)
    
    exponential_decay(Xj_endpoints[period_by_timesteps], 
                      target_by_timesteps, 
                      const_by_timesteps, 
                      countup_by_timesteps)
  }
  
  out
  
}

simulate_outcome = function(timesteps, Bt, sum_Xj, Z_init = 0, eta, sampling_timestep) {
  
  # Draw process noise for simulation steps
  Z_noise = rnorm(length(timesteps)-1, sd = sd_p * sqrt(sampling_timestep))
  
  B = Bt(timesteps)
  Zt = numeric(length(timesteps))
  Zt[1] = Z_init
  
  for (i in 1:(length(timesteps) - 1)) {
    Zt[i+1] = exponential_decay(Zt[i], B[i] + sum_Xj[i], 
                                eta, sampling_timestep) + Z_noise[i]
  }
  
  approxfun(timesteps, Zt, method = 'linear')
}