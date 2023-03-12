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