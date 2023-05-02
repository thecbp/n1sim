#' Simulate data from an N-of-1 trial with washout
#'
#' This function simulates data from an N-of-1 trial, which is a single-subject 
#' clinical trial in which the subject serves as their own control. The function
#'  generates a dataset that includes time, period, treatment, baseline, 
#'  treatment effects, simulated outcome without error, and simulated outcome 
#'  with error.
#'
#' @param n_trts The number of treatments being tested.
#' @param n_blocks The number of blocks, where each block consists of a sequence of treatments.
#' @param period_length The length of each treatment period.
#' @param washout_length The length of the washout period between treatments, during which no treatment is administered.
#' @param sampling_timestep The time interval at which observations are sampled.
#' @param sd_b Standard deviation of the baseline process.
#' @param sd_p Standard deviation of the process noise.
#' @param sd_o Standard deviation of the observation noise.
#' @param Ejs A vector of treatment effect sizes for each treatment (Ej).
#' @param gammas A vector of decay constants for the decay of treatment effects when the treatment is not being applied (gamma).
#' @param taus A vector of decay constants for the decay of treatment effects when the treatment is being applied (tau).
#' @param eta Time constant for the decay of the effect of the total treatment on the outcome (Zt).
#' @param phi AR constant for the autocorrelation of the underyling outcome process
#' @param custom.Bt Custom baseline function.
#' @param Xj_inits Custom initial values for treatment effect functions.
#' @param Z_init Custom initial value for the outcome process (Zt).
#'
#' @return a list with three main components: "data", "functions", and "params"
#' - \code{data}: A tibble (data frame) containing the simulated data for the N-of-1 trial with columns:
#'    - code{time}: Observation time points.
#'    - code{period}: Period number corresponding to each time point.
#'    - code{treatment}: Treatment administered at each time point.
#'    - code{Bt}: Baseline values at each time point.
#'    - code{sum_Xjt}: Total treatment effect at each time point.
#'    - code{Zt}: Simulated outcome without error at each time point.
#'    - code{Yt}: Simulated outcome with error at each time point.
#'    - code{Tj}: Binary indicators for each treatment at each time point (one column per treatment).
#'    - code{Xj}: Treatment effect for each treatment at each time point (one column per treatment).
#' - \code{functions} A list of functions used during the simulation:
#'    - \code{Bt}: Baseline function
#'    - \code{Tt}: Function mapping time points to treatments.
#'    - \code{Xjt}: List of treatment effect functions for each treatment.
#'    - \code{Zt}: Simulated outcome function without error.
#'  -\code{params}: A list of parameter settings:
#'    - \code{trial}: List of trial parameters, including number of treatments (n_trts), number of blocks (n_blocks), period length (period_length), and washout length (washout_length).
#'    - \code{process}: List of process parameters, such as time steps (timesteps), sampling time step (sampling_timestep), and standard deviations for baseline (sd_b), process noise (sd_p), and observation noise (sd_o).
#'    - \code{treatment}: List of treatment effect parameters, including treatment effect sizes (Ejs), decay constants for treatment effects when treatment is not being applied (gammas), decay constants for treatment effects when treatment is being applied (taus), and the time constant for the decay of the effect of the total treatment on the outcome (eta).
#'
#' @importFrom magrittr %>%
#'
#' @export simulate


simulate = function(n_trts, n_blocks, period_length, washout_length,
                    sampling_timestep, sd_b, sd_p, sd_o,
                    Ejs, gammas, taus, eta, phi,
                    custom.Bt = NULL, Xj_inits = NULL, Z_init = NULL) {
  
  # TO-DO:
  # incorporate an interaction effects?
  
  # Preventing checks to these variables in documentation
  if (FALSE) {
    . <- NULL
    time <- NULL
  }
  
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
      incorporate_washout(., n_trts)
    switches = find_treatment_switches(n_blocks * n_trts, period_length, washout_length)
  }
  
  # Construct necessary functions for simulations
  Bt = simulate_baseline_function(timesteps, sd_b, sampling_timestep, phi)
  
  # Treatment number by timestep
  Tt = stats::approxfun(switches, c(order, order[length(order)]), method = 'constant')
  
  # Construct the effect functions for each treatment
  Xjt = lapply(1:n_trts, function(j) { 
    Tjt = as.numeric(Tt(timesteps) == j)
    simulate_effect_function(j, order, Tjt, switches, Xj_init = 0, Ejs[j], taus[j], gammas[j] )
  })
  
  # Get the total treatment effect across time
  sum_Xj = lapply(Xjt, function(f) { f(timesteps) }) %>% Reduce("+", .)
  
  Zt = simulate_outcome(timesteps, Bt, sum_Xj, Z_init = 0, eta, sampling_timestep, sd_p)
  Yt = Zt(timesteps) + stats::rnorm(length(timesteps), 0, sd_o)
  
  data = tibble::tibble(
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
incorporate_washout = function(trt_order, n_trts) {
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

simulate_baseline_function = function(time_vec, sd_b, dt, phi) {
  
  # The original Percha simulation incorporates drift into the baseline function
  # but I will replace this with an AR1 process
  
  
  # Generate outcome based on if we need serial correlation or not
  if (phi == 0) {
    proc = stats::rnorm(length(time_vec) - 1, 0, sd_b * sqrt(dt))
  } else {
    proc = stats::arima.sim(model = list(ar = phi), sd = sd_b * sqrt(dt), n = length(time_vec) - 1)
  }
  
  baseline = c(time_vec[1], proc)
  stats::approxfun(time_vec, baseline, method = 'linear')
}

simulate_treatment_function = function(n_trts, n_blocks, trt_changes) {
  # counterbalanced + randomized treatment order
  # add n_blocks at end to enable function interpolation
  block_perms = combinat::permn(1:n_trts)
  block_vec = c(sample(block_perms, size = n_blocks, replace = T) %>% unlist(), n_blocks)
  
  list(
    block_vec = block_vec,
    Tt = stats::approxfun(trt_changes, block_vec, method = 'constant')
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

simulate_outcome = function(timesteps, Bt, sum_Xj, Z_init = 0, eta, sampling_timestep, sd_p) {
  
  # Draw process noise for simulation steps
  Z_noise = stats::rnorm(length(timesteps)-1, sd = sd_p * sqrt(sampling_timestep))
  
  B = Bt(timesteps)
  Zt = numeric(length(timesteps))
  Zt[1] = Z_init
  
  for (i in 1:(length(timesteps) - 1)) {
    Zt[i+1] = exponential_decay(Zt[i], B[i] + sum_Xj[i], 
                                eta, sampling_timestep) + Z_noise[i]
  }
  
  stats::approxfun(timesteps, Zt, method = 'linear')
}