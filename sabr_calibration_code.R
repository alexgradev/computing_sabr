# ========================================================================
# ATM SABR implied volatility approximation (special-case Hagan formula)
#
# Purpose: Calculate the model implied volatility of a ATM option using 
#          the special-case formula
#
# Inputs:
#  --- F_0:   forward price/rate (> 0); numeric scalar or vector
#  --- H:     time to expiry in years (> 0); numeric scalar
#  --- alpha: SABR volatility level (> 0); numeric scalar
#  --- beta:  SABR elasticity (typically in [0, 1]); numeric scalar
#  --- rho:   correlation parameter with abs(rho) < 1; numeric scalar
#  --- v:     vol-of-vol (>= 0); numeric scalar
#
# Output:
#  --- vol_atm: implied Black volatility at-the-money; numeric vector

sabr_atm_vol <- function(F_0, H, alpha, beta, rho, v) {
  alpha / (F_0^(1 - beta)) *
    (1 + H * (
      ((1 - beta)^2 / 24) * (alpha^2 / (F_0^(2 - 2*beta))) +
        ((rho * beta * v) / 4) * (alpha / (F_0^(1 - beta))) +
        ((2 - 3 * rho^2) / 24) * (v^2)
    ) 
    )
}

  
# ========================================================================
# SABR implied volatility approximation (general Hagan formula)
#
# Purpose: Computes SABR model implied volatility
#
# Inputs:
#  --- F_0:   forward price/rate (F_0 > 0); numeric scalar or vector
#  --- K:     strike (K > 0); numeric scalar or vector
#  --- H:     time to expiry in years (H > 0); numeric scalar
#  --- alpha: SABR volatility level (alpha > 0); numeric scalar
#  --- beta:  SABR elasticity (typically in [0, 1]); numeric scalar
#  --- rho:   correlation parameter with abs(rho) < 1; numeric scalar
#  --- v:     vol-of-vol (v >= 0); numeric scalar
#
# Output:
#  --- vol: implied Black volatility for each (F_0, K) pair; numeric vector
#
# Notes:
#  --- Uses the function sabr_atm_vol(F_0, H, alpha, beta, rho, v) 
#      for the ATM case

sabr_vol <- function(F_0, K, H, alpha, beta, rho, v) {
  
  if (any(F_0 <= 0) || any(K <= 0)) stop("F_0 and K must be > 0 
                                     (Black/lognormal SABR).")
  if (alpha <= 0)                 stop("alpha must be > 0")
  if (H <= 0)                     stop("H must be > 0")
  if (abs(rho) >= 1)              stop("|rho| must be < 1")
  if (v < 0)                     stop("v must be >= 0")
  
  # Define variables for repeatedly appearing calculations
  one_m_b <- 1 - beta
  logFK   <- log(F_0 / K)    # log-moneyness
  abslog  <- abs(logFK)
  FK_beta <- (F_0 * K)^(one_m_b / 2)
  
  # Denominator correction - adjusts for log(F_0/K) terms
  log2 <- logFK^2
  log4 <- log2^2
  denom <- FK_beta * (1 + (one_m_b^2/24) * log2 + (one_m_b^4/1920) * log4)
  
  # Hagan variables z and x(z)
  z  <- (v / alpha) * FK_beta * logFK
  xz <- log((sqrt(1 - 2*rho*z + z^2) + z - rho) / (1 - rho))
  
  # z/x(z) has a removable singularity at z=0 ATM
  # Use the limit z/x(z) -> 1 for numerical stability
  z_over_x <- ifelse(abs(z) < 1e-08, 1, z / xz)
  
  # Time correction term - accounts for expiry horizon H
  term1 <- (one_m_b^2/24) * (alpha^2 / ((F_0*K)^(one_m_b)))
  term2 <- (rho*beta*v/4) * (alpha / FK_beta)
  term3 <- ((2 - 3*rho^2) * v^2 / 24)
  
  vol <- (alpha / denom) * z_over_x * (1 + (term1 + term2 + term3) * H)
  
  # Near ATM situation -> use explicit ATM formula when log(F_0/K) is very small
  idx_atm <- abslog < 1e-10
  if (any(idx_atm)) vol[idx_atm] <- sabr_atm_vol(F_0, H, alpha, beta, rho, v)
  
  vol
}

    

# ========================================================================
# Numerical rootfinding of alpha from ATM via the cubic (Helper)
#
# Purpose: If we rearrange the ATM approximation equation to solve for alpha, 
#          we obtain a polynomial equation in alpha. In this implementation the 
#          rearrangement is written as a cubic:
#               a3 * alpha^3 + a2 * alpha^2 + a1 * alpha + a0 = 0
#
#          This function constructs the cubic coefficients (a0..a3) implied by 
#          the ATM approximation and returns *all* roots of the polynomial 
#          using polyroot().
#
# Inputs:
#  --- rho:       SABR correlation (abs(rho) < 1)
#  --- v:         vol-of-vol (>= 0)
#  --- F_0:       forward price/rate (> 0)
#  --- H:         time to expiry in years (> 0)
#  --- beta:      elasticity (typically in [0, 1])
#  --- sigma_atm: observed market ATM implied Black volatility (> 0)
#
# Output:
#  --- roots: complex vector of length 3 containing the cubic roots for alpha
#
# Notes:
#  --- This is a numerical helper, real_alpha_roots() is used to pick alpha
#  --- polyroot() can return complex roots even when real solutions exist
#
# ========================================================================

alpha_roots <- function(rho, v, F_0, H, beta, sigma_atm) {
  
  # a3*alpha^3 + a2*alpha^2 + a1*alpha + a0 = 0
  a3 <- (1 - beta)^2 * H / (24 * F_0^(2 - 2*beta))
  a2 <- rho * beta * v * H / (4 * F_0^(1 - beta))
  a1 <- 1 + (2 - 3 * rho^2) * v^2 * H / 24
  a0 <- -sigma_atm * F_0^(1 - beta)
  
  # Find the roots
  polyroot(c(a0, a1, a2, a3))
}

# ========================================================================
# Numerical finding of real roots from ATM via the cubic 
#
# Purpose: Calls alpha_roots(...) to compute all cubic roots for alpha.
#          Then it filters those roots to keep only roots that are (numerically)
#          real (imaginary part close to zero) and positive. Finally, it returns
#          the smallest positive real root
#
# Inputs:
#  --- rho:       SABR correlation (|rho| < 1); controls skew/asymmetry
#  --- v:         vol-of-vol (>= 0); controls smile curvature
#  --- F_0:       forward price/rate (> 0)
#  --- H:         time to expiry in years (> 0)
#  --- beta:      elasticity (typically in [0, 1])
#  --- sigma_atm: observed market ATM implied Black volatility (> 0)
#  --- imag_tol:  tolerance for treating complex roots as real (default 1e-6)
#
# Output:
#  --- alpha_hat: calibrated alpha implied by ATM, chosen as the smallest 
#                 positive real root of the cubic
#
# Notes:
#  --- Stops with an error if no positive real root is found (handled later
#      in the full calibration routine via tryCatch(...) block)
#  --- Numerical root solvers typically return tiny imaginary parts due to 
#      floating point error (e.g. 0.5000000 + 1e-12 i). We treat such roots as
#      real if the imaginary part is below a tolerance parameter imag_tol
#
# ========================================================================

real_alpha_roots <- function(rho, v, F_0, H, beta, sigma_atm, imag_tol = 1e-6) {
  # Find all roots
  roots <- alpha_roots(rho, v, F_0, H, beta, sigma_atm)
  
  # Take only real positive roots
  real_roots          <- Re(roots)[abs(Im(roots)) < imag_tol]
  positive_real_roots <- real_roots[real_roots > 0]
  
  if (length(positive_real_roots) == 0) stop("No positive real root for alpha.")
  
  # Return the smallest root
  min(positive_real_roots)
}


# ========================================================================
# Year fraction between two dates
#
# Purpose: Convert two calendar dates into a time-to-expiry expressed in years.
#          This is needed in SABR because the model formulas use a continuous 
#          time horizon H in years
#
# Inputs:
#  --- SettleDate:   valuation/settlement date (Date / POSIXct)
#  --- ExerciseDate: option expiry/exercise date (Date / POSIXct)
#
# Output:
#  --- H: time difference in years (numeric)
#
# Notes:
#  --- Uses 365.242199 days per year (mean tropical year) as a constant scaling
#
# ========================================================================

yearfrac <- function(SettleDate, ExerciseDate) {
  as.numeric(difftime(ExerciseDate, SettleDate, units = "days")) / 365.242199
}
    

# ========================================================================
# Calibration of SABR parameters (rho, v) for ONE expiry by using implied
# alpha from ATM
#
# Purpose: Calibrate the SABR smile for a single maturity H by:
#          (1) identifying the ATM quote (closest strike to F_0)
#          (2) solving for alpha from the ATM volatility (via real_alpha_roots),
#              given candidate (rho*, v*) and fixed beta
#          (3) fitting (hat_rho, hat_nu) by minimizing the sum of squared errors 
#              between market implied vols and SABR model vols across strikes
#          (4) recomputing the final alpha implied by ATM using the fitted 
#              (hat_rho, hat_nu).
#
# Inputs:
#  --- SettleDate:          settlement date (used in yearfrac)
#  --- ExerciseDate:        option expiry/exercise date (used in yearfrac)
#  --- F_0:                 forward price/rate at SettleDate (> 0)
#  --- K_vec:               vector of strikes (> 0), same length as vol_vec
#  --- vol_vec:             vector of market implied Black vols, same length as
#                           K_vec
#  --- beta:                fixed SABR elasticity (usually in [0, 1], often 0.5)
#  --- start:               initial guess for optimizer parameters (rho, v)
#  --- lower:               lower bounds for (rho, v) in optimization
#  --- upper:               upper bounds for (rho, v) in optimization
#  --- include_atm_in_fit:  if FALSE (default), sets the ATM residual to 0 so ATM
#                           is not overweighted (alpha is already pinned to ATM)
#
# Output (list):
#  --- Alpha:    final alpha implied from ATM using calibrated (hat_rho, hat_nu)
#  --- Beta:     beta used (fixed)
#  --- Rho:      calibrated rho (smile skew/asymmetry)
#  --- Nu:       calibrated v (smile curvature / vol-of-vol)
#  --- SSE:      final sum of squared errors across strikes
#  --- conv:     optim() convergence code; 0 typically indicates success
#  --- idx_atm:  index of the strike treated as ATM (closest to F_0)
#
# Notes:
#  --- Uses yearfrac(SettleDate, ExerciseDate) to compute H in years
#  --- Uses optim(..., method="L-BFGS-B") to fit (rho, v) with bounds
#  --- Uses penalty residuals (1e3) to keep optimizer away from invalid regions
#      (abs(rho) close to 1, v < 0 or when alpha cannot be solved robustly)
#  --- Requires sabr_vol(), real_alpha_roots() and their dependencies.
#
# ========================================================================

calibrate_sabr_one_expiry <- function(SettleDate, ExerciseDate, F_0, K_vec, vol_vec,
                                      beta = 0.5,
                                      start = c(0, 0.5),
                                      lower = c(-0.999, 1e-8),
                                      upper = c( 0.999, 10.0),
                                      include_atm_in_fit = FALSE) {
  
  # Identify which strike is ATM (closest to F_0)
  idx_atm <- which.min(abs(K_vec - F_0))
  sigma_atm <- vol_vec[idx_atm]
  
  # Convert dates to year fraction
  H <- yearfrac(SettleDate, ExerciseDate)
  
  # Residuals vector (MarketVols - ModelVols)
  residual_vec <- function(X) {
    rho <- X[1]
    v  <- X[2]
    
    # Boundary penalties - keep optimizer away from invalid regions
    if (abs(rho) > 0.9999 || v < 0) return(rep(1e3, length(K_vec)))
    
    # Robust calculation of the real positive roots of alpha
    alpha <- tryCatch(
      real_alpha_roots(rho, v, F_0, H, beta, sigma_atm),
      error = function(e) NA_real_
    )
    # Penalize when alpha cannot be solved robustly
    if (!is.finite(alpha) || alpha <= 0) return(rep(1e3, length(K_vec)))
    
    # Compute model vols across strikes using Hagan SABR
    model <- sabr_vol(F_0, K_vec, H, alpha, beta, rho, v)
    
    # Residuals (vector)
    res <- vol_vec - model
    
    # (OPTIONAL) Remove ATM from fitting because alpha is forced to match ATM to
    # prevent giving the ATM extra influence relative to the other strikes 
    if (!include_atm_in_fit) res[idx_atm] <- 0
    
    res
  }
  
  # Least-squares objective
  # minimize sum(residuals^2)
  sse_obj <- function(X) {
    r <- residual_vec(X)
    sum(r^2)
  }
  
  # Optimization over (rho, v) with bounds
  fit <- optim(
    par    = start,
    fn     = sse_obj,
    method = "L-BFGS-B",
    lower  = lower,
    upper  = upper
  )
  
  # Extract calibrated parameters
  rho_hat <- fit$par[1]
  nu_hat  <- fit$par[2]
  
  # Final alpha implied from ATM using calibrated (rho_hat, nu_hat)
  alpha_hat <- real_alpha_roots(rho_hat, nu_hat, F_0, H, beta, sigma_atm)
  
  list(
    Alpha   = alpha_hat,
    Beta    = beta,
    Rho     = rho_hat,
    Nu      = nu_hat,
    SSE     = fit$value,
    conv    = fit$convergence,
    idx_atm = idx_atm
  )
}

# ========================================================================
# Calibration of a SABR volatility SURFACE across multiple expiries
#
# Purpose: Build a SABR parameter surface by calibrating SABR independently for
#          each expiry (each column of the market smile data). For every expiry:
#            (1) pick the forward F_0 for that maturity
#            (2) take the vector of strikes K and market vols for that maturity
#            (3) call calibrate_sabr_one_expiry(...) to fit (rho, v) and imply
#                alpha from the ATM quote (given fixed beta)
#          The results are then collected into a parameter table (one row per 
#          expiry)
#
# Inputs:
#  --- SettleDate:           settlement date (used in yearfrac)
#  --- ExerciseDates:        vector of expiries (length = n_exp)
#  --- MarketStrikes:        matrix of strikes with dimensions [n_strikes x n_exp]
#                            where column j corresponds to expiry ExerciseDates[j]
#  --- MarketVolatilities:   matrix of market implied Black vols with the same
#                            dimensions as MarketStrikes
#  --- CurrentForwardValues: vector of forwards (length = n_exp), one forward per
#                            expiry
#  --- beta:                 fixed SABR elasticity (typically in [0, 1], often 0.5)
#  --- start:                initial guess for (rho, v) for each expiry calibration
#  --- include_atm_in_fit:   passed through to calibrate_sabr_one_expiry(); if FALSE
#                            (default), ATM residual is set to 0 to avoid overweighting
#
# Output (list):
#  --- params:  data.frame with one row per expiry containing:
#               Alpha, Beta, Rho, Nu, SSE, conv
#  --- details: list of length n_exp with the full calibration output for each expiry
#               (including idx_atm and any other diagnostics returned)
#
# Notes:
#  --- This is a per-expiry calibration, i.e. no smoothing across maturities
#  --- Assumes MarketStrikes and MarketVolatilities are aligned column-wise to
#      ExerciseDates and CurrentForwardValues
#  --- Requires calibrate_sabr_one_expiry() and all underlying helpers
#
# ========================================================================

calibrate_sabr_surface <- function(SettleDate, ExerciseDates,
                                   MarketStrikes, MarketVolatilities,
                                   CurrentForwardValues,
                                   beta = 0.5,
                                   start = c(0, 0.5),
                                   include_atm_in_fit = FALSE) {
  
  n_exp <- length(ExerciseDates)
  res_list <- vector("list", n_exp)
  
  # Iterate over each expiry to obtain the surface
  for (j in 1:n_exp) {
    res_list[[j]] <- calibrate_sabr_one_expiry(
      SettleDate         = SettleDate,
      ExerciseDate       = ExerciseDates[j],
      F_0                  = as.numeric(CurrentForwardValues[j]),
      K_vec              = as.numeric(MarketStrikes[, j]),
      vol_vec            = as.numeric(MarketVolatilities[, j]),
      beta               = beta,
      start              = start,
      include_atm_in_fit = include_atm_in_fit
    )
  }
  
  # Pack results into a convenient table
  out <- data.frame(
    Alpha = sapply(res_list, function(x) x[["Alpha"]]),
    Beta  = sapply(res_list, function(x) x[["Beta"]]),
    Rho   = sapply(res_list, function(x) x[["Rho"]]),
    Nu    = sapply(res_list, function(x) x[["Nu"]]),
    SSE   = sapply(res_list, function(x) x[["SSE"]]),
    conv  = sapply(res_list, function(x) x[["conv"]])
  )
  
  list(params = out, details = res_list)
}

# ========================================================================
# Plot market volatility smile vs SABR fitted smile for one expiry 
#
# Purpose: Visualize how well the calibrated SABR parameters match the observed
#          market implied Black volatilities for a single expiry index j.
#          The function takes the j-th maturity slice, computes the SABR fitted
#          volatilities using the calibrated parameters and overlays the market
#          smile (points) and SABR fitted smile (line)
#
# Inputs:
#  --- j:                    integer index of the expiry slice to plot
#  --- results:              calibration output from calibrate_sabr_surface()
#                            (must contain results$details[[j]] 
#                            with Alpha/Beta/Rho/Nu)
#  --- Settle:               settlement/valuation date (used in yearfrac)
#  --- ExerciseDates:        vector of expiries (length = n_exp)
#  --- MarketStrikes:        matrix [n_strikes x n_exp] of strikes 
#                            (column j = expiry j)
#  --- MarketVolatilities:   matrix [n_strikes x n_exp] of market implied Black 
#                            vols aligned with MarketStrikes
#  --- CurrentForwardValues: vector of forwards (length = n_exp), one per expiry
#
# Output:
#  --- A base R plot
# ========================================================================

plot_smile <- function(j, results, Settle, ExerciseDates, MarketStrikes,
                       MarketVolatilities, CurrentForwardValues) {
  
  pars    <- results$details[[j]]
  
  F_0     <- as.numeric(CurrentForwardValues[j])
  K       <- as.numeric(MarketStrikes[, j])
  vol_mkt <- as.numeric(MarketVolatilities[, j])
  H       <- yearfrac(Settle, ExerciseDates[j])
  
  vol_fit <- sabr_vol(F_0, K, H, pars$Alpha, pars$Beta, pars$Rho, pars$Nu)
  
  plot(K, vol_mkt, pch = 16, xlab = "Strike", ylab = "Black vol",
       main = paste0(Maturities[j]," into ",Tenor, " Volatility Smile"))
  
  lines(K, vol_fit, lwd = 2)
  
  legend("topright", legend = c("Market", "SABR fit"),
         pch = c(16, NA), lwd = c(NA, 2), bty = "n")
}

# ============================================================================

# ========================================================================
# Plot SABR volatility surface
#
# Purpose:
#   Plot the calibrated SABR implied volatility surface as a 3D surface using
#   base R's persp(). The surface is evaluated on a common strike grid across
#   all expiries, while maturities are represented by year fractions H
#
# What it does:
#   1) Builds H_vec = yearfrac(SettleDate, ExerciseDates) (x-axis)
#   2) Builds a common strike grid K_grid (y-axis)
#   3) For each expiry j, pulls SABR params (Alpha, Beta, Rho, Nu) from fit$details[[j]]
#      and computes vols on K_grid using sabr_vol(F0, K_grid, H_j, ...)
#   4) Assembles the z-matrix for persp() with orientation:
#        x = H_vec  (length nH)
#        y = K_grid (length nK)
#        z must be a matrix with nrow(z) = nH and ncol(z) = nK
#      Therefore we compute vol_mat as [nK x nH] and pass z = t(vol_mat).
#
# Inputs:
#  --- results:              calibration output from calibrate_sabr_surface() 
#                            (must contain results$details list)
#  --- SettleDate:           Date object used to compute maturities in years
#  --- ExerciseDates:        Date object used to compute maturities in years
#  --- MarketStrikes:        strike matrix [n_strikes x n_exp] 
#                            (used only to set grid range)
#  --- CurrentForwardValues: vector of forwards per expiry (length n_exp)
#  --- strike_grid:          optional numeric vector; 
#                            if NULL, a grid is created automatically
#  --- n_strike:             number of strike points if strike_grid is NULL
#  --- theta, phi, expand, ticktype: passed to persp() for viewing angle/scale
#
# Output:
#   - Draws a 3D surface plot
# ========================================================================


plot_surface <- function(results,
                         SettleDate, ExerciseDates,
                         MarketStrikes, CurrentForwardValues,
                         strike_grid = NULL,
                         n_strike = 60,
                         theta = 45, phi = 25, expand = 0.7,
                         ticktype = "detailed") {
  
  n_exp <- length(ExerciseDates)
  if (n_exp != ncol(MarketStrikes))     
                    stop("ExerciseDates length must match ncol(MarketStrikes).")
  
  if (n_exp != length(CurrentForwardValues)) 
                    stop("CurrentForwardValues length must match ExerciseDates.")
  
  if (is.null(results$details) || length(results$details) != n_exp) 
                    stop("results must be output of calibrate_sabr_surface().")
  
  # x-axis: maturities in years
  H_vec <- yearfrac(SettleDate, ExerciseDates)
  
  # y-axis: strike grid 
  if (is.null(strike_grid)) {
    K_min <- min(MarketStrikes, na.rm = TRUE)
    K_max <- max(MarketStrikes, na.rm = TRUE)
    K_grid <- seq(K_min, K_max, length.out = n_strike)
  } else {
    K_grid <- as.numeric(strike_grid)
  }
  
  # Compute vols on grid: vol_mat_KH is [K x H]
  vol_mat_KH <- matrix(NA_real_, nrow = length(K_grid), ncol = n_exp)
  for (j in 1:n_exp) {
    pars <- results$details[[j]]
    F0   <- as.numeric(CurrentForwardValues[j])
    Hj   <- H_vec[j]
    
    vol_mat_KH[, j] <- sabr_vol(F0, K_grid, Hj, pars$Alpha, pars$Beta, pars$Rho, pars$Nu)
  }
  
  # persp needs z as [H x K] because x=H_vec and y=K_grid
  z_HK <- t(vol_mat_KH)
  
  persp(x = H_vec, y = K_grid, z = z_HK,
        xlab = "Maturity (years)", ylab = "Strike", zlab = "Black vol",
        main = "SABR Implied Volatility Surface", 
        theta = theta, phi = phi, expand = expand, ticktype = ticktype)
}

# ========================================================================
# ========================================================================
# ========================================================================

SettleDate     <- as.Date('16-Jan-2026', format = "%d-%b-%Y")

ExerciseDatesChar <- c('16-Mar-2026','16-Jan-2027','16-Jan-2028','16-Jan-2029',
                       '16-Jan-2030','16-Jan-2031','16-Jan-2033','16-Jan-2036')
ExerciseDates     <- as.Date(ExerciseDatesChar, format = "%d-%b-%Y")

YearsToExercise <- yearfrac(SettleDate, ExerciseDates)
Maturities      <- c("3M", "1Y", "2Y", "3Y", "4Y", "5Y", "7Y", "10Y")
Tenor           <- "10Y"
MarketVolatilities <- matrix(c(
  # 3M,   1Y,   2Y,   3Y,   4Y,   5Y,   7Y,   10Y
  57.6, 53.7, 49.4, 45.6, 44.1, 41.1, 35.2, 32.0,
  46.6, 46.9, 44.8, 41.6, 39.8, 37.4, 33.4, 31.0,
  35.9, 39.3, 39.6, 37.9, 37.2, 34.7, 30.5, 28.9,
  34.1, 36.5, 37.8, 36.6, 35.0, 31.9, 28.1, 26.6, # ATM 
  41.0, 41.3, 39.5, 37.8, 36.0, 32.6, 29.0, 26.0,
  45.8, 43.4, 41.9, 39.2, 36.9, 33.2, 29.6, 26.3,
  50.3, 46.9, 44.0, 40.0, 37.5, 33.8, 30.2, 27.3
), nrow = 7, byrow = TRUE) / 100  

MarketStrikes <- matrix(c(
  # 3M,   1Y,   2Y,   3Y,   4Y,   5Y,   7Y,   10Y
  1.00, 1.25, 1.68, 2.00, 2.26, 2.41, 2.58, 2.62,
  1.50, 1.75, 2.18, 2.50, 2.76, 2.91, 3.08, 3.12,
  2.00, 2.25, 2.68, 3.00, 3.26, 3.41, 3.58, 3.62,
  2.50, 2.75, 3.18, 3.50, 3.76, 3.91, 4.08, 4.12, # ATM 
  3.00, 3.25, 3.68, 4.00, 4.26, 4.41, 4.58, 4.62,
  3.50, 3.75, 4.18, 4.50, 4.76, 4.91, 5.08, 5.12,
  4.00, 4.25, 4.68, 5.00, 5.26, 5.41, 5.58, 5.62
), nrow = 7, byrow = TRUE) / 100

# Current Forward Rates and ATM Vols across the maturities
CurrentForwardValues <- MarketStrikes[4,]
ATMVolatilities      <- MarketVolatilities[4,]

results_calibration <- calibrate_sabr_surface(
  SettleDate = SettleDate,
  ExerciseDates = ExerciseDates,
  MarketStrikes = MarketStrikes,
  MarketVolatilities   = MarketVolatilities,
  CurrentForwardValues = CurrentForwardValues,
  beta  = 0.5,
  start = c(0, 0.5),
  include_atm_in_fit = FALSE)

plot_smile(j = 1,
           results_calibration, 
           SettleDate, ExerciseDates,
           MarketStrikes, 
           MarketVolatilities, 
           CurrentForwardValues)

# plot_surface(results_calibration,
#              SettleDate,
#              ExerciseDates,
#              MarketStrikes,
#              CurrentForwardValues,
#              n_strike = 80)






